from submine.io.transcode import detect_format, FMT_GSPAN, FMT_LG, FMT_EDGELIST


def test_detect_format_gspan_data_variants():
    assert detect_format("foo.data") == FMT_GSPAN
    assert detect_format("foo.data.x") == FMT_GSPAN
    assert detect_format("foo.data.2") == FMT_GSPAN
    assert detect_format("foo.DATA.3") == FMT_GSPAN


def test_detect_format_lg_and_edgelist():
    assert detect_format("graph.lg") == FMT_LG
    assert detect_format("graph.edgelist") == FMT_EDGELIST
    assert detect_format("graph.txt") == FMT_EDGELIST
