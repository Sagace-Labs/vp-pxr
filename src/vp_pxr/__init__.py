"""Pregnane X receptor activation from a SMILES string.

A ligand in the receptor's pocket induces CYP3A4 and the transporters that
accompany it, which redirects the clearance of anything co-administered and can
raise the flux through a bioactivating route.

    from vp_pxr import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])   # -> DataFrame[pxr_agonist, pxr_cytotox]
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from vp_core.registry import Version, VersionedPathway
from vp_pxr.target import CYTOTOX, TARGET, TARGETS, Endpoint, all_names
from vp_pxr.target import get as get_target

__version__ = "1.0.0"

PATHWAY = "pxr"
VERSIONS_DIR = Path(__file__).resolve().parent / "versions"


def _predict_values(model: Any, smiles: list[str], version: Version) -> np.ndarray:
    """Columns for ``version``, in the order its signature declares them."""
    from rdkit import Chem, RDLogger

    from vp_core import fingerprints, xgb

    RDLogger.DisableLog("rdApp.*")
    X = fingerprints.featurize(smiles, str(version.features))
    values = np.column_stack(
        [xgb.predict_proba(model[name], X) for name in version.output_names]
    ).astype(np.float32)

    # An unparseable input is a declared NaN.
    unparseable = [Chem.MolFromSmiles(s) is None for s in smiles]
    values[np.asarray(unparseable)] = np.nan
    return values


_pathway = VersionedPathway(PATHWAY, VERSIONS_DIR, predict_fn=_predict_values)

__all__ = [
    "CYTOTOX",
    "PATHWAY",
    "TARGET",
    "TARGETS",
    "VERSIONS_DIR",
    "Endpoint",
    "Version",
    "__version__",
    "all_names",
    "current_version",
    "get",
    "get_target",
    "predict",
    "signature",
    "versions",
]


def predict(smiles: list[str], *, version: str | None = None) -> pd.DataFrame:
    """Score each SMILES with ``version`` (default: newest)."""
    return _pathway.predict(smiles, version=version)


def versions() -> list[str]:
    """Released version names, oldest first."""
    return _pathway.versions()


def current_version() -> str:
    """The newest released version."""
    return _pathway.current()


def get(version: str | None = None) -> Version:
    """Load a version and its validated manifest."""
    return _pathway.get(version)


def signature(version: str | None = None) -> dict[str, Any]:
    """The output contract of a version."""
    return get(version).manifest["signature"]
