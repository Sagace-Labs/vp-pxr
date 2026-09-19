"""The output contract.

``train`` stamps this into each new version's manifest.

Bump ``SIGNATURE_VERSION`` only for a **breaking** change — a removed or
renamed column, a changed dtype, or a changed meaning. Adding a new column is
also breaking for anyone consuming positionally.
"""

from __future__ import annotations

from typing import Any

__all__ = [
    "INPUTS",
    "OUTPUTS",
    "PRIMARY",
    "SIGNATURE_VERSION",
    "as_manifest_table",
    "column_names",
]

SIGNATURE_VERSION = 1

INPUTS: list[str] = ["smiles"]

OUTPUTS: list[dict[str, Any]] = [
    {
        "name": "pxr_agonist",
        "dtype": "float32",
        "range": [0.0, 1.0],
        "semantics": (
            "P(activates the human pregnane X receptor reporter in the Tox21 qHTS "
            "agonist screen)"
        ),
        "missing": "NaN when RDKit cannot parse the input SMILES",
    },
    {
        "name": "pxr_cytotox",
        "dtype": "float32",
        "range": [0.0, 1.0],
        "semantics": (
            "P(reduces viability in the counter-screen over the same library). A "
            "compound scoring high on both moved the reporter in a cell that was "
            "also dying"
        ),
        "missing": "NaN when RDKit cannot parse the input SMILES",
    },
]

#: The endpoint the version's headline metrics describe.
PRIMARY = OUTPUTS[0]["name"]


def column_names() -> list[str]:
    return [o["name"] for o in OUTPUTS]


def as_manifest_table() -> dict[str, Any]:
    """The ``[signature]`` table as it is written into a manifest."""
    return {
        "version": SIGNATURE_VERSION,
        "inputs": list(INPUTS),
        "outputs": [dict(o) for o in OUTPUTS],
    }
