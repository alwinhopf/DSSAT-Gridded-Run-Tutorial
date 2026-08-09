# One-point DSSAT model smoke fixture

This fixture is intentionally small and is committed under the workspace test
fixture exception in `AGENTS.md` section 5. It contains only the four essential
point inputs; the verified DSSAT installation supplies standard/genotype data
through its adjacent `DSSATPRO` file.

- Location: 43.102691 N, 92.937244 W (Iowa)
- Weather: NASA POWER daily data for 1992, written in DSSAT `.WTH` format
- Soil: SSURGO-derived profile formatted for DSSAT
- Experiment: `UFGA9201.MZX`, adapted to point/weather/soil ID `00000001`
- Expected DSSAT 4.8 result: `HWAM = 6292 kg/ha`, `CWAM = 15971 kg/ha`

Run explicitly on a machine with DSSAT installed:

```bash
DSSAT_EXE=/verified/path/to/dscsm048 \
  python -m pytest tests/test_dssat_model_e2e.py -v
Rscript --vanilla tests/test_dssat_model_e2e.R
```

The tests skip when no executable or adjacent `DSSATPRO.V48` /
`DSSATPRO.L48` is available. Hosted provider CI keeps model validation separate
and reports that DSSAT was unavailable rather than pretending it executed.
