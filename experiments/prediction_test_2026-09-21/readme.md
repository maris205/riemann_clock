# Prospective Fe II relative-frequency prediction test

This directory contains a **design, not new observations or a registered protocol**. It asks whether a universal mean relative-frequency response can predict unseen absorbers after its amplitude has been estimated from a separate development sample. Such a mean drift is an additional hypothesis beyond the original Riemann-accessibility uncertainty model.

The primary descriptor is `D = v(Fe II 2382) - v(Fe II 2374)`. For the explicitly fixed cosmology, reference time and inverse-log-square schedule, illustrative redshifts 1.0, 1.5 and 2.0 give:

```text
D(1.5) = 0.457563556 D(1.0) + 0.542436444 D(2.0).
```

A constant also obeys this relation. An independently repeated nonzero endpoint contrast, calibrated nuisance control and sufficient precision are necessary. No amplitude, sign or atomic response has been derived from the Riemann model. All 100 m/s examples are hypothetical. Powers 1, 2 and 3 remain nearly indistinguishable at practical precision.

- [Full Chinese protocol](prediction_protocol_cn.md)
- [Instrument feasibility](reports/instrument_feasibility_cn.md)
- [Independent prediction-logic review](reports/prediction_logic_review.md)
- [Machine-readable calculations](results/design.json)
- [Design figure](figures/prediction_design.pdf)

Reproduce locally:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/prediction_design.py
```

The existing science fits and manuscripts remain unchanged. Actual targets, exposure times, independent calibration constraints, the prediction amplitude, and a public registration remain to be established before a confirmatory observation campaign.
