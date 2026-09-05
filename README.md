# ILuRF SPICE files builder

## Requirements

- Linux

- Wget, GNU Fortran, GCC

    Assuming you are on Debian or similar distribution:
  
    ```
    sudo apt-get install wget gfortran gcc
    ```

- Calceph 5.0.0 or later
  
    If you are on Debian 13+ or similar distribution:
    
    ```
    sudo apt-get install libcalceph-dev
    ```

    For older or non-Debian distributions, please see the [installation instructions](https://calceph.imcce.fr/installation/linux).
  
## Download ephemeris

```
./download-epm.sh &&
./download-inpop.sh &&
cd de430 && ./download-de430.sh
```

### Convert DE430 to the binary DE format

(Assuming you are still in the `de430` directory)

```
./convert-de430.sh
cd ..
```

## Build

```
./build.sh
```

## Run

```
./ilurf-spice
```

The output files are `ilurf2026.bsp` and `ilurf2026.bpc`.

## License

The entire contents of this repository is public domain.
