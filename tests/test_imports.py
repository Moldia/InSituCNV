def test_package_import_smoke():
    import insitucnv

    assert hasattr(insitucnv, "run_insitucnv")
