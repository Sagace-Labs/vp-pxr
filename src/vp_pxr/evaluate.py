"""Measure a version under an evaluation protocol.

    python -m vp_pxr.evaluate --version v1

Writes ``versions/<version>/metrics.json`` and regenerates ``CARD.md``. The
protocol named in the version's manifest supplies the split, the fold sizes,
the seed set and the metric list. ``--protocol`` measures the same weights
under another registered protocol; ``metrics.json`` keys every run by
protocol id.

Each seed refits on its own training fold.

Every output is scored on the same folds, restricted to the compounds its own
endpoint labels. The first declared output carries the headline metrics.

A seed that yields an empty fold raises rather than being skipped.
"""

from __future__ import annotations

import argparse
import sys
from datetime import date

import numpy as np

from vp_core import fingerprints
from vp_pxr import contract as pxr_contract
from vp_pxr import data as pxr_data
from vp_pxr import model as pxr_model
from vp_pxr.train import OUTPUT_LABELS

__all__ = ["evaluate_version", "main"]


def _score_output(X, table, column, smiles, protocol) -> tuple[dict, list[dict]]:
    """Per-seed metrics for one endpoint, on the folds its compounds fall in."""
    from vp_core import metrics as metrics_mod

    rows = pxr_data.labelled(table, column)
    keep = set(rows.tolist())
    y = table[column].to_numpy()

    per_seed: list[dict[str, float]] = []
    folds: list[dict[str, int]] = []
    for seed in protocol.seeds:
        train_idx, val_idx, test_idx = protocol.split_indices(smiles, seed)
        train = np.array([i for i in train_idx if i in keep])
        val = np.array([i for i in val_idx if i in keep])
        test = np.array([i for i in test_idx if i in keep])
        if not len(train) or not len(val) or not len(test):
            raise ValueError(
                f"seed {seed} leaves {column!r} with an empty fold — this endpoint "
                "cannot be measured under this protocol"
            )

        fitted = pxr_model.fit(
            X[train], y[train].astype(int), X[val], y[val].astype(int), seed=seed
        )
        from vp_core import xgb

        proba = xgb.predict_proba(fitted, X[test])
        scored = metrics_mod.binary_metrics(y[test].astype(int), proba)
        per_seed.append(scored)
        folds.append(
            {
                "seed": int(seed),
                "train": len(train),
                "val": len(val),
                "test": len(test),
            }
        )
        print(
            f"  {column} seed {seed}: test AUC {scored['auc_roc']:.4f} "
            f"(n_test={len(test)})",
            file=sys.stderr,
        )

    return (
        {
            "n": {
                "total": len(rows),
                "train": folds[0]["train"],
                "val": folds[0]["val"],
                "test": folds[0]["test"],
            },
            "folds": folds,
            "test": metrics_mod.aggregate(per_seed, protocol.metrics),
        },
        folds,
    )


def evaluate_version(
    version: str,
    *,
    protocol_id: str | None = None,
    use_example: bool = False,
    write: bool = True,
) -> dict:
    """Run the version's protocol and return the metrics record."""
    import vp_pxr
    from vp_core import card, dataset, metrics_store, protocols
    from vp_core import manifest as manifest_mod

    resolved = vp_pxr.get(version)
    protocol_id = protocol_id or resolved.protocol_id
    protocol = protocols.get(protocol_id)
    labels = manifest_mod.dataset_labels(resolved.manifest)

    table = pxr_data.example() if use_example else pxr_data.load()
    smiles = table["smiles"].tolist()
    X = fingerprints.featurize(smiles, pxr_model.FEATURES)

    scored = {
        output: _score_output(X, table, OUTPUT_LABELS[output], smiles, protocol)[0]
        for output in resolved.output_names
    }
    primary = scored[pxr_contract.PRIMARY]

    record = {
        "pathway": "pxr",
        "version": resolved.name,
        "protocol_id": protocol_id,
        "dataset_sha256": dataset.dataset_hash(table, labels=labels),
        "dataset": "example fixture" if use_example else "full",
        "evaluated": date.today().isoformat(),
        "n": primary["n"],
        "folds": primary["folds"],
        "test": primary["test"],
        "additional_outputs": [
            {"name": name, **scored[name]}
            for name in resolved.output_names
            if name != pxr_contract.PRIMARY
        ],
    }

    if write:
        if use_example:
            raise ValueError(
                "refusing to record fixture metrics as a released result — "
                "--example is for smoke-checking the harness only"
            )
        path = metrics_store.write(resolved.directory, record)
        card.write_card(resolved.directory)
        auc = record["test"]["auc_roc"]
        print(
            f"wrote {path}\n"
            f"  {protocol_id}: AUC {auc['mean']:.4f} +/- {auc['std']:.4f} "
            f"over {len(protocol.seeds)} seeds",
            file=sys.stderr,
        )
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m vp_pxr.evaluate")
    parser.add_argument("--version", default=None, help="default: the newest version")
    parser.add_argument(
        "--protocol",
        default=None,
        help="default: the protocol the version's manifest declares",
    )
    parser.add_argument(
        "--example",
        action="store_true",
        help="run on the committed fixture without writing (smoke check only)",
    )
    args = parser.parse_args(argv)

    import vp_pxr

    version = args.version or vp_pxr.current_version()
    record = evaluate_version(
        version,
        protocol_id=args.protocol,
        use_example=args.example,
        write=not args.example,
    )
    if args.example:
        auc = record["test"]["auc_roc"]["mean"]
        print(f"fixture smoke AUC {auc if np.isfinite(auc) else float('nan'):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
