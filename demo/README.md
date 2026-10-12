# Demo

Run `scorch demo` from the repository folder after installing SCORCH
(see [docs/installation.md](../docs/installation.md)). It writes a synthetic
example to `demo/run/` and runs the whole workflow on it in a few seconds.
[docs/demo.md](../docs/demo.md) explains the output.

Files in this folder:

- `demo_overview.png`: the figure produced with the recorded seed (synthetic
  data, not observations)
- `expected_results.json`: the events the demo must find; checked by
  `tests/test_demo.py`
- `fit_power_law_demo.R`: the optional R step of `scorch demo --with-r`
