#!/bin/sh
for n in `seq 1550 100 2550`; do
    wget -nc --no-check-certificate https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de430/ascp$n.430;
done;
wget -nc --no-check-certificate https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de430/header.430_229;
wget -nc --no-check-certificate https://ssd.jpl.nasa.gov/ftp/eph/planets/fortran/asc2eph.f

