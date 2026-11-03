import random
from pathlib import Path

import pytest

from spot.config import Gepard as GepardOptions
from spot.dotplot import GepardPlotter
from spot.tools import find_gepard
from spot.validate import Region

gepard = find_gepard()
pytestmark = pytest.mark.skipif(gepard is None, reason="Gepard not installed")


def test_many_dotplots_in_few_jvms(tmp_path: Path) -> None:
    assert gepard is not None
    rng = random.Random(3)
    jobs = []
    for i in range(5):
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        a = Region(f"a{i}", "chr1", 1, 3000, seq, 1000, 1000)
        b = Region(f"b{i}", "chr1", 1, 3000, seq[::-1], 1000, 1000)
        jobs.append((tmp_path / f"run/og{i}/dotplots/a{i}__b{i}.png", a, b))
    assert GepardPlotter(gepard, GepardOptions(), threads=2).plot(jobs) == 5
    for out, _, _ in jobs:
        assert out.is_file() and out.stat().st_size > 0
