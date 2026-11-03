from pathlib import Path

import pytest

from spot import config
from spot.errors import ConfigError


def test_defaults(config_file: Path) -> None:
    cfg = config.load(config_file)
    assert cfg.profile == "default"
    assert cfg.minimap2.x == "asm10"
    assert cfg.anchors.types == ["SYN", "INV"]
    assert cfg.validation.enabled


def test_relative_paths_follow_the_config_file(config_file: Path) -> None:
    cfg = config.load(config_file)
    assert cfg.orthofinder.root == config_file.parent / "orthofinder"
    assert cfg.species["AA"].genome == config_file.parent / "genomes" / "a.fasta"
    assert cfg.species["BB"].genome == Path("/abs/b.fasta")
    assert cfg.outdir == config_file.parent / "spot_out"


def test_profile_then_file_then_overrides(config_file: Path) -> None:
    cfg = config.load(config_file, profile="relaxed")
    assert cfg.anchors.types == ["SYN", "INV", "TRANS", "INVTR", "DUP", "INVDP"]

    config_file.write_text(
        config_file.read_text() + "profile: relaxed\nvalidation: {min_frac: 0}\nanchors:\n  types: [SYN]\n"
    )
    assert config.load(config_file).anchors.types == ["SYN"]

    overrides = {"anchors": {"types": "syn,inv,syn"}, "minimap2": {"x": "asm20"}}
    cfg = config.load(config_file, overrides=overrides)
    assert cfg.anchors.types == ["SYN", "INV"]
    assert cfg.minimap2.x == "asm20"
    assert cfg.validation.min_frac == 0


def test_assignment() -> None:
    assert config.assignment("validation.min_frac=0.3") == {"validation": {"min_frac": 0.3}}
    assert config.assignment("anchors.types=[SYN, TRANS]") == {"anchors": {"types": ["SYN", "TRANS"]}}
    with pytest.raises(ConfigError):
        config.assignment("validation.min_frac")


def test_native_tool_options_pass_through(config_file: Path) -> None:
    cfg = config.load(config_file, overrides={"syri": {"invgaplen": 500}, "blastn": {"evalue": "1e-5"}})
    assert cfg.syri.model_dump() == {"nosnp": True, "invgaplen": 500}
    assert cfg.blastn.evalue == 1e-5


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"anchors": {"types": ["SYN", "FOO"]}}, "unknown SyRI class 'FOO'"),
        ({"minimap2": {"x": "asm7"}}, "minimap2.x"),
        ({"query": "AA"}, "must be different"),
        ({"reference": "ZZ"}, "reference 'ZZ'"),
        ({"validation": {"min_frac": 2}}, "validation.min_frac"),
        ({"anchors": {"gap": 1}}, "anchors.gap"),
        ({"species": {"AA": {"match": "(["}}}, "invalid regex"),
    ],
)
def test_invalid(config_file: Path, overrides: dict, message: str) -> None:
    with pytest.raises(ConfigError, match=message.replace("(", r"\(").replace("[", r"\[")):
        config.load(config_file, overrides=overrides)


def test_species_count(config_file: Path) -> None:
    extra = {"species": {"CC": {"label": "C", "match": "^CC_"}}}
    with pytest.raises(ConfigError, match="exactly two species"):
        config.load(config_file, overrides=extra)


def test_unknown_profile(config_file: Path) -> None:
    with pytest.raises(ConfigError, match="available: default, relaxed"):
        config.load(config_file, profile="strictest")


def test_profile_file(config_file: Path) -> None:
    (config_file.parent / "mine.yaml").write_text("description: x\nanchors: {types: [SYN, DUP]}\n")
    assert config.load(config_file, profile="mine.yaml").anchors.types == ["SYN", "DUP"]
    (config_file.parent / "bad.yaml").write_text("outdir: elsewhere\n")
    with pytest.raises(ConfigError, match="profiles only set"):
        config.load(config_file, profile="bad.yaml")


def test_profiles_skip_appledouble(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "profiles").mkdir()
    (tmp_path / "profiles" / "default.yaml").write_text("description: d\n")
    (tmp_path / "profiles" / "._default.yaml").write_bytes(b"\x00\x05\x16\x07\xb0")
    monkeypatch.setattr(config.resources, "files", lambda package: tmp_path)
    assert list(config.builtin_profiles()) == ["default"]


def test_builtin_profiles_are_valid() -> None:
    for name in config.builtin_profiles():
        settings = config.profile_settings(name)
        assert settings["minimap2"]["x"] == "asm10"
