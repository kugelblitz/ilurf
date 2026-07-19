/*
 * This program builds a lunar ephemeris (geocentric position and libration)
 * averaged from DE430, EPM2021, and INPOP21a using predefined weights. The
 * result is stored in SPICE files (BSP for position, binary PCK for libration)
 * and can be used as a realization of the International Lunar Reference Frame.
 *
 * Dependency: Calceph library (available in a Linux distribution near you).
 * Author: Dmitry Pavlov, 2026.
 * Public domain.
 */

#include <stdio.h>
#include <calceph.h>
#include <stdlib.h>
#include <math.h>

/* Macros and global variables are awesome */

double weights[3] = 
  {
    0.167596331684500, // DE
    0.451454072137127, // EPM
    0.380949596178372  // INPOP
  };

// Number of Chebyshev coefficients stored per record
// (same for orbit and libration)
#define NUM_COEF  15

// Size of grid on which the Chebyshev approximation is made in one record
// (all knots belong to (-1, 1))
#define NUM_KNOTS 16 

// Precomputed grid knots and values of Chebyshev polynomials in them
static double chebyshev_knots[NUM_KNOTS];
static double chebyshev_values_for_knots[NUM_KNOTS][NUM_COEF];

/* Helper functions for Chebyshev approximation */

void chebyshev_values (double x, double T[NUM_COEF])
{
  T[0] = 1;
  T[1] = x;

  for (int k = 2; k < NUM_COEF; k++)
    T[k] = 2 * x * T[k - 1] - T[k - 2];
}

void calc_chebyshev_coefficients (double f[NUM_KNOTS], double coef[NUM_COEF])
{
  // c_0 = 1/n * \sum f(x_i)
  coef[0] = 0.0;
  for (int i = 0; i < NUM_KNOTS; i++)
    coef[0] += f[i];
  coef[0] /= NUM_KNOTS;
  
  // c_k = 2/n * \sum f(x_i) * T_k(x_i)
  for (int k = 1; k < NUM_COEF; k++)
  {
    coef[k] = 0.0;
    for (int i = 0; i < NUM_KNOTS; i++)
      coef[k] += f[i] * chebyshev_values_for_knots[i][k];
    coef[k] *= 2.0 / NUM_KNOTS;
  }
}


int main (void)
{
  /* Open input files */
  
  t_calcephbin *eph[3];
  const char* filenames_epm[] = {"epm2021.bsp", "moonlibr_epm2021.bpc", "moonlibr_epm2021.tf"};
  const char* filenames_inpop[] = {"inpop21a_TDB_m1000_p1000_spice.bsp", "inpop21a_TDB_m1000_p1000_spice.bpc",
                                   "inpop21a_TDB_m1000_p1000_spice.tf"};

  eph[0] = calceph_open("de430/JPLEPH");
  eph[1] = calceph_open_array(3, filenames_epm);
  eph[2] = calceph_open_array(3, filenames_inpop);
  
  if (!eph[0] || !eph[1] || !eph[2])
  {
    fprintf(stderr, "error opening planetary ephemeris file for reading\n");
    return -1;
  }

  double jd_from = 2440587.5; // 1970-01-01
  double jd_to   = 2470587.5; // 2052-02-20
  double step_orbit     = 2.0; // orbit is stored by 2-day subintervals
  double step_libration = 5.0; // libration is stored by 5-day intervals
  int    nrecords_orbit     = (int)((jd_to - jd_from) / step_orbit);
  int    nrecords_libration = (int)((jd_to - jd_from) / step_libration);
  double jd_cur;

  /* Open output files */
  t_writephbin *weph_orbit     = writeph_spk_create("ilurf2026.bsp", "ILuRF2026_Orbit", 0);
  t_writephbin *weph_libration = writeph_pck_create("ilurf2026.bpc", "ILuRF2026_Libration", 0);

  if (!weph_orbit || !weph_libration)
  {
    fprintf(stderr, "error opening planetary ephemeris file for writing\n");
    return -1;
  }
  
  /* Precalculate Chebyshev stuff */
  for (int i = 0; i < NUM_KNOTS; i++)
  {
    chebyshev_knots[i] = cos(M_PI * (i + 0.5) / NUM_KNOTS);
    chebyshev_values(chebyshev_knots[i], chebyshev_values_for_knots[i]);
  }

  // All Chebyshev coefficients that are to be written into the orbit file
  // (a lesser part of the array will be reused for the libration file)
  double *polynomials = (double *)malloc(sizeof(double) * nrecords_orbit * NUM_COEF * 3);
    
  /* Read orbits, calculate Chebyshev coefficients */
  for (int irecord = 0; irecord < nrecords_orbit; irecord++)
  {
    double p_knots[3][NUM_KNOTS];

    for (int i = 0; i < NUM_KNOTS; i++)
    {
      double pv[3][6];
      
      // Map knot from [-1, 1] to the timespan of the current record
      double jd_knot = jd_from + step_orbit * (irecord + 0.5 * (chebyshev_knots[i] + 1));

      // Get positions and velocities from the three ephemerides
      for (int j = 0; j < 3; j++)
      {
        if (!calceph_compute_unit(eph[j], jd_knot, 0, 10, 3, CALCEPH_UNIT_KM + CALCEPH_UNIT_SEC, pv[j]))
        {
          fprintf(stderr, "error accessing lunar position\n");
          return -1;
        }
      }

      // Weighted sum
      for (int j = 0; j < 3; j++)
        p_knots[j][i] = pv[0][j] * weights[0] + pv[1][j] * weights[1] + pv[2][j] * weights[2];
    }

    for (int j = 0; j < 3; j++)
      calc_chebyshev_coefficients(p_knots[j], polynomials + NUM_COEF * (3 * irecord + j));
  }

  // Save orbit
  if (!writeph_spk2_seq_write(weph_orbit, NAIFID_MOON, NAIFID_EARTH, 1,
                              jd_from, 0.0, jd_to, 0.0,
                              step_orbit, polynomials, nrecords_orbit,
                              NUM_COEF - 1, "Moon"))
  {
    fprintf(stderr, "error writing lunar orbit\n");
    return -1;
  }

  writeph_close(weph_orbit);
                           
  printf("Stored %d records of orbit\n", nrecords_orbit);
  
  /* Read librations, calculate Chebyshev coefficients */
  
  for (int irecord = 0; irecord < nrecords_libration; irecord++)
  {
    int i;
    double a_knots[3][NUM_KNOTS];

    for (int i = 0; i < NUM_KNOTS; i++)
    {
      double libr[3][6];
    
      // Map knot from [-1, 1] to the timespan of the current record
      double jd_knot = jd_from + step_libration * (irecord + 0.5 * (chebyshev_knots[i] + 1));

      // Get Euler angles and their rates from the three ephemerides
      for (int j = 0; j < 3; j++)
      {
        if (!calceph_orient_unit(eph[j], jd_knot, 0, NAIFID_MOON,
                                 CALCEPH_UNIT_RAD + CALCEPH_UNIT_SEC + CALCEPH_USE_NAIFID, libr[j]))
        {
          fprintf(stderr, "error accessing lunar libration\n");
          return -1;
        }
      }

      // INPOP has different zero for longitude libration;
      // adjust to match others in the weighted sum
      libr[2][2] += 2 * M_PI * 408;
      
      // Weighted sum
      for (int j = 0; j < 3; j++)
        a_knots[j][i] = libr[0][j] * weights[0] + libr[1][j] * weights[1] + libr[2][j] * weights[2];
    }

    for (int j = 0; j < 3; j++)
      calc_chebyshev_coefficients(a_knots[j], polynomials + NUM_COEF * (3 * irecord + j));
  }

  // Save orbit
  if (!writeph_pck2_seq_write(weph_libration, NAIFID_MOON, 1,
                              jd_from, 0.0, jd_to, 0.0,
                              step_libration, polynomials, nrecords_libration,
                              NUM_COEF - 1, "moon"))
  {
    fprintf(stderr, "error writing lunar libration\n");
    return -1;
  }

  writeph_close(weph_libration);
  
  printf("Stored %d subintevals of libration\n", nrecords_libration);

  calceph_close(eph[0]);
  calceph_close(eph[1]);
  calceph_close(eph[2]);
  free(polynomials);
}
