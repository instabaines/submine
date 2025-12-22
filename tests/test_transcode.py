from __future__ import annotations

from pathlib import Path

from submine.io.transcode import transcode_path, FMT_GSPAN


def data_path(name: str) -> Path:
    return Path(__file__).parent / "data" / name


def test_transcode_edgelist_to_gspan(tmp_path: Path):
    src = data_path("sample.edgelist")
    dst = tmp_path / "out.data"

    transcode_path(src, dst, dst_fmt=FMT_GSPAN)

    lines = dst.read_text(encoding="utf-8").strip().splitlines()
    assert lines == [
        "t # 0",
        "v 0 2",
        "v 1 2",
        "v 2 2",
        "e 0 1 2",
        "e 1 2 2",
        "e 2 0 2",
        "t # -1",
    ]
