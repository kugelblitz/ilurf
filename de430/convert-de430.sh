wget --no-check-certificate https://ssd.jpl.nasa.gov/ftp/eph/planets/fortran/asc2eph.f &&
sed -i '/PARAMETER ( NRECL = 4 )/s/^C//' asc2eph.f &&
gfortran asc2eph.f -o asc2eph
cat header.430_229 ascp*.430 | ./asc2eph
