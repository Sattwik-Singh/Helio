"""Scoring detections against a reference catalogue.

A detection matches a catalogue event if their times differ by no more than
`tolerance_s`; each catalogue event can be claimed by at most one detection
(closest in time wins). From the matches:

    precision = matched detections / all detections
    recall    = matched catalogue events / all catalogue events

Both are needed: a detector that fires constantly has perfect recall, and one
that fires once, correctly, has perfect precision.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MatchResult:
    pairs: list[tuple[int, int]]  # (detection index, catalogue index)
    false_positives: list[int]
    misses: list[int]
    precision: float
    recall: float
    f1: float
    median_abs_time_error_s: float


def match_events(detected_times, catalogue_times, tolerance_s: float = 300.0) -> MatchResult:
    det = np.asarray(detected_times, dtype=float)
    cat = np.asarray(catalogue_times, dtype=float)
    candidates = [(abs(d - c), i, j) for i, d in enumerate(det) for j, c in enumerate(cat) if abs(d - c) <= tolerance_s]
    candidates.sort()
    used_d, used_c, pairs = set(), set(), []
    for _, i, j in candidates:
        if i in used_d or j in used_c:
            continue
        used_d.add(i)
        used_c.add(j)
        pairs.append((i, j))
    fp = [i for i in range(len(det)) if i not in used_d]
    miss = [j for j in range(len(cat)) if j not in used_c]
    precision = len(pairs) / len(det) if len(det) else float("nan")
    recall = len(pairs) / len(cat) if len(cat) else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if len(pairs) else 0.0
    err = float(np.median([abs(det[i] - cat[j]) for i, j in pairs])) if pairs else float("nan")
    return MatchResult(pairs=sorted(pairs), false_positives=fp, misses=miss, precision=precision, recall=recall, f1=f1,
                       median_abs_time_error_s=err)
