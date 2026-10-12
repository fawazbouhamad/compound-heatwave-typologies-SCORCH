# Input folder

Put the processed temperature file here as `input/tmax_processed.nc`.

For the study, download `SCORCH_data.zip` from the Zenodo dataset record
(version 1.0.0, DOI [10.5281/zenodo.21717874](https://doi.org/10.5281/zenodo.21717874)),
check its MD5 (`80772df9410306cbc419191d0dcb824c`), extract it and copy
`scorch_data/tmax_processed.nc` into this folder. The expected SHA-256 of the
NetCDF file is `fc080880b6f2065fb5eec9fe9dce71abdfdc44a72d1ae91e50abdd278fadaa90`.

Check the file with:

```bash
scorch check-input input/tmax_processed.nc --strict
```

Everything else in this folder is ignored by Git. To try SCORCH without the
dataset, run `scorch demo`.
