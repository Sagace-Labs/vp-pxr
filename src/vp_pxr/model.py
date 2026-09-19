"""The PXR model: binary XGBoost on Morgan bits, MACCS keys and descriptors.

Hyperparameters and feature choice.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from vp_core import fingerprints, xgb

__all__ = ["FEATURES", "HYPERPARAMS", "fit", "predict"]

FEATURES = "combo3"

HYPERPARAMS: dict[str, Any] = {
    "n_estimators": 2000,
    "learning_rate": 0.05,
    "max_depth": 6,
    "min_child_weight": 1.0,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "reg_lambda": 2.0,
    "reg_alpha": 0.05,
    "gamma": 0.1,
    "early_stopping_rounds": 40,
}


def fit(X_train, y_train, X_val, y_val, *, seed: int = 0) -> Any:
    """Fit one endpoint's binary model. The validation fold stops boosting."""
    return xgb.fit_binary(
        X_train,
        y_train,
        X_val,
        y_val,
        params=dict(HYPERPARAMS),
        seed=seed,
    )


def predict(model: Any, smiles: list[str]) -> np.ndarray:
    """Positive-class probability for arbitrary SMILES."""
    return xgb.predict_proba(model, fingerprints.featurize(smiles, FEATURES))
