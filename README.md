# vp-pxr

Predicts activation of the pregnane X receptor from a SMILES string — a
molecular initiating event for drug-induced liver injury, where a ligand in the
receptor's pocket induces CYP3A4 and the transporters that accompany it.

## Install

    pip install vp-pxr

## Use

    from vp_pxr import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])   # -> DataFrame[pxr_agonist, pxr_cytotox]

Returns one row per input and one column per declared output. `pxr_agonist` is
the probability of activating the receptor's reporter; `pxr_cytotox` is the
probability of reducing viability in the counter-screen over the same library.
A compound scoring high on both moved the reporter in a cell that was also
dying. Unparseable SMILES come back as NaN. Pin a version with `predict(smiles,
version="v1")`; list what is available with `versions()`.

## Current version

**v1**, signature 1, measured under protocol `scaffold-balanced-5seed@1`.
The full record — metrics per output and per seed, dataset hash,
environment — is in
[`src/vp_pxr/versions/v1/CARD.md`](src/vp_pxr/versions/v1/CARD.md).

## Data

PubChem BioAssay AID 1346982, the Tox21 qHTS screen for agonists of the human
pregnane X receptor, and AID 1346977, its cell-viability counter-screen over the
same library. Both are reduced to one row per compound labelled by the majority
call across its assay records, retrieved 2026-09-06 and redistributed here as a
United States government work in the public domain. Rebuild and check for
upstream drift with `python -m vp_pxr.data fetch --verify`; see
[`data/README.md`](data/README.md) for the expected layout.

## Retrain

    python -m vp_pxr.train --version v2 --reason "why this version exists"
    python -m vp_pxr.evaluate --version v2

`train` fits one deployment model per output on the whole dataset and writes a
new version directory; `evaluate` refits per seed under the protocol and
records what those held-out models scored. Reproducibility is to the recorded
dataset hash and environment, which can change.

## Licence

Code is Apache-2.0 ([`LICENSE`](LICENSE)). The bundled dataset is in the public
domain ([`LICENSE-DATA`](LICENSE-DATA)).

## Cite

See [`CITATION.cff`](CITATION.cff).
