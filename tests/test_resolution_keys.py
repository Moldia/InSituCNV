"""The CNV Leiden cluster key convention is shared across modules.

``cluster_cnv_resolutions`` and ``run_insitucnv`` must agree on the obs key so
``primary_resolution`` points at a column that exists.
"""

from __future__ import annotations

from insitucnv.tl.cnv import cnv_leiden_key


def test_cnv_leiden_key_format():
    assert cnv_leiden_key(0.1) == "cnv_leiden_res0.1"
    assert cnv_leiden_key(0.25) == "cnv_leiden_res0.25"
    assert cnv_leiden_key(1.0) == "cnv_leiden_res1"
