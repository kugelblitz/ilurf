#!/bin/sh
for F in epm2021.bsp epm2021.tpc moonlibr_epm2021.bpc moonlibr_epm2021.tf; do   
    wget --no-check-certificate https://ftp.iaaras.ru/pub/epm/EPM2021/SPICE/$F;
done

# add a line for Calceph to read the lunar ephemeris without problems
echo "      OBJECT_MOON_FRAME = 'MOON_PA_EPM2021'" >> moonlibr_epm2021.tf

