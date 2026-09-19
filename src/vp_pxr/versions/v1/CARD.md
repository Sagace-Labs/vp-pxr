# pxr v1

_Generated from `manifest.toml` and `metrics.json`. Do not edit._

Released 2026-09-19 · signature 1

**Why this version.** first release: binary XGBoost on Morgan, MACCS and RDKit descriptors over the Tox21 pregnane-X-receptor screen

## Outputs

| column | dtype | range | meaning |
|---|---|---|---|
| `pxr_agonist` | float32 | 0.0–1.0 | P(activates the human pregnane X receptor reporter in the Tox21 qHTS agonist screen) |
| `pxr_cytotox` | float32 | 0.0–1.0 | P(reduces viability in the counter-screen over the same library). A compound scoring high on both moved the reporter in a cell that was also dying |

Missing values: NaN when RDKit cannot parse the input SMILES

## Performance

Protocol `scaffold-balanced-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed and each group placed in the fold it overfills least, so a group larger than a fold's capacity settles in train instead of starving that fold. Same fold fractions, seeds and metrics as scaffold-shuffle-5seed@1; only the packing differs. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-09-19 on n_train=5274, n_val=483, n_test=725.

### `pxr_agonist`

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.835 | 0.017 | 0.811, 0.845, 0.819, 0.844, 0.856 |
| auprc | 0.707 | 0.035 | 0.651, 0.692, 0.711, 0.757, 0.726 |
| mcc | 0.501 | 0.025 | 0.454, 0.518, 0.495, 0.512, 0.526 |
| brier | 0.171 | 0.013 | 0.190, 0.166, 0.179, 0.168, 0.153 |

### `pxr_cytotox`

Measured on n_train=4918, n_val=425, n_test=624.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.869 | 0.016 | 0.873, 0.889, 0.876, 0.841, 0.863 |
| auprc | 0.557 | 0.024 | 0.570, 0.571, 0.586, 0.519, 0.540 |
| mcc | 0.463 | 0.048 | 0.499, 0.524, 0.456, 0.383, 0.453 |
| brier | 0.101 | 0.014 | 0.113, 0.090, 0.084, 0.120, 0.098 |

> Comparable only with metrics carrying the same protocol id.

## Data

Tox21 pregnane X receptor qHTS — PubChem BioAssay AID 1346982 (Human pregnane X receptor (PXR) small molecule agonists, qHTS assay) and AID 1346977 (Human pregnane X receptor (PXR) small molecule agonists, qHTS cell viability counter screen), rows called Active or Inactive, one row per compound labelled by majority call across its assay records. Retrieved 2026-09-06, licensed public-domain, redistributed here.

`6482` compounds, positive rate `0.244`, table SHA-256 `eea5ba76176f73b4…`

Regenerate and check for upstream drift with `python -m vp_pxr.data fetch --verify`.

## Model

xgboost-binary on `combo3` features. Shipped weights: one model per output, each on every compound its endpoint labels minus a 10% scaffold carve used for early stopping

`weights.joblib` SHA-256 `006712c723fda23b…`

## Provenance

Environment: python 3.11.11, rdkit 2026.03.5, xgboost 3.2.0.

Reproducibility is to this dataset hash and this environment, not bit-exact: the sources are live endpoints and RDKit descriptor values move between releases.
