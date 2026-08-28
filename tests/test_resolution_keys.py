"""Regression: the cluster-key convention must be shared across modules.

``find_optimal_clustering`` used to write ``cnv_leiden_<r>`` while
``select_best_resolution`` / ``run_insitucnv`` looked up ``cnv_leiden_res<r>``,
so ``select_resolution_by_metrics=True`` pointed at a non-existent obs column.
"""

from __future__ import annotations

import inspect

from insitucnv import analysis
from insitucnv.tl.cnv import cnv_leiden_key


def test_cnv_leiden_key_format():
    assert cnv_leiden_key(0.1) == "cnv_leiden_res0.1"
    assert cnv_leiden_key(0.25) == "cnv_leiden_res0.25"
    assert cnv_leiden_key(1.0) == "cnv_leiden_res1"


def test_find_optimal_clustering_uses_canonical_key():
    src = inspect.getsource(analysis.find_optimal_clustering)
    assert "cnv_leiden_key(" in src
    assert 'f"cnv_leiden_{' not in src
