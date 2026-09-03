# ILuRF SPICE files builder

## Requirements

- Linux

- Wget, GNU Fortran, Calceph, GCC

    Assuming you are on Debian 13+ or similar distribution:
    
    ```
    sudo apt-get install wget gfortran gcc libcalceph-dev
    ```

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
