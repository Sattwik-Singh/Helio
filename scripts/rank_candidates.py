"""A learned ranking of candidates: logistic regression on six audit features.

    python scripts/rank_candidates.py results/benchmark_<train>_detections.csv results/benchmark_<test>_detections.csv [...]

Trained on the first file, evaluated on each of the others without refitting.
The label is "is in the catalogue", which is imperfect: a real shock that the
catalogue lacks counts as a negative. The model is deliberately tiny (seven
numbers) because the training year has only a few dozen events.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize

FEATURES = ["log_b_ratio", "log_n_ratio", "mach_fast_low", "sharpness", "log_normal_sigma", "triggers"]


def features(t: pd.DataFrame) -> np.ndarray:
    f = pd.DataFrame({
        "log_b_ratio": np.log(t.b_ratio), "log_n_ratio": np.log(t.n_ratio),
        "mach_fast_low": t.mach_fast_low.clip(0, 4), "sharpness": t.sharpness.clip(-0.5, 1.5),
        "log_normal_sigma": np.log(t.normal_sigma_deg.clip(0.5, 90)), "triggers": t.triggers.clip(0, 30),
    })
    return f[FEATURES].fillna(f[FEATURES].median()).to_numpy()


def fit(x, y, l2=1.0):
    mu, sd = x.mean(0), x.std(0) + 1e-9
    z = np.column_stack([np.ones(len(x)), (x - mu) / sd])

    def loss(w):
        s = z @ w
        return np.mean(np.logaddexp(0, s) - y * s) + l2 * np.sum(w[1:] ** 2) / len(y)

    w = minimize(loss, np.zeros(z.shape[1]), method="BFGS").x
    return lambda xx: 1 / (1 + np.exp(-np.column_stack([np.ones(len(xx)), (xx - mu) / sd]) @ w)), w


def auc(score, y):
    pos, neg = score[y == 1], score[y == 0]
    return float(np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg]))


def report(name, score, y, n_catalogue, thresholds):
    print(f"\n{name}: {int(y.sum())} catalogued among {len(y)} candidates ({n_catalogue} in catalogue), AUC {auc(score, y):.3f}")
    print("  threshold  kept  matched  false+  precision  recall")
    for th in thresholds:
        keep = score >= th
        tp, fp = int((keep & (y == 1)).sum()), int((keep & (y == 0)).sum())
        print(f"  {th:>9.2f}  {keep.sum():>4}  {tp:>7}  {fp:>6}  {tp / max(tp + fp, 1):>9.3f}  {tp / n_catalogue:>6.3f}")


def main(paths):
    train = pd.read_csv(paths[0])
    y = train.matched.to_numpy().astype(float)
    model, w = fit(features(train), y)
    print("weights on standardised features:", dict(zip(["bias"] + FEATURES, np.round(w, 2))))
    thresholds = (0.2, 0.3, 0.5, 0.7)
    import json

    n_catalogue = lambda p: json.load(open(p.replace("_detections.csv", "_summary.json")))["catalogue_shocks"]  # noqa: E731
    report(f"TRAIN {paths[0]}", model(features(train)), y, n_catalogue(paths[0]), thresholds)
    for p in paths[1:]:
        t = pd.read_csv(p)
        n_cat = n_catalogue(p)
        report(f"HELD OUT {p}", model(features(t)), t.matched.to_numpy().astype(float), n_cat, thresholds)


if __name__ == "__main__":
    main(sys.argv[1:])
