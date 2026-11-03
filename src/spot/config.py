from __future__ import annotations

import os
import re
from collections.abc import Mapping
from importlib import resources
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator, model_validator

from .errors import ConfigError

SYRI_CLASSES = ("SYN", "INV", "TRANS", "INVTR", "DUP", "INVDP")
HEADER_REGEX = r"(?P<chr>[^\s:]+):(?P<start>\d+)-(?P<end>\d+)"


class _Section(BaseModel):
    model_config = ConfigDict(extra="forbid")


class _ToolOptions(BaseModel):
    model_config = ConfigDict(extra="allow")


def _f(default: Any, text: str, **kw: Any) -> Any:
    """A field with a description (shown as help in the app)."""
    if callable(default):
        return Field(default_factory=default, description=text, **kw)
    return Field(default, description=text, **kw)


class Minimap2(_ToolOptions):
    x: Literal["asm5", "asm10", "asm20"] = _f(
        "asm10", "minimap2 preset: asm5 for ~0.1 % divergence, asm10 ~1 %, asm20 ~5 %"
    )


class Syri(_ToolOptions):
    nosnp: bool = _f(True, "skip SNP and indel calling (not needed by SPOT, much faster)")


class Blastn(_ToolOptions):
    evalue: float = _f(1e-10, "maximum E-value of a flank HSP", gt=0)
    perc_identity: float = _f(80.0, "minimum identity (%) of a flank HSP", ge=0, le=100)


class Gepard(_ToolOptions):
    matrix: str = _f("edna.mat", "substitution matrix (a name in Gepard's matrices folder, or a path)")
    word: int = _f(10, "word length of the suffix-array dotplot", ge=1)


class Anchors(_Section):
    types: list[str] = _f(lambda: ["SYN", "INV"], "SyRI classes whose alignments are used as anchors")
    same_chromosome: bool = _f(True, "only anchors joining homologous chromosomes")
    homology: Literal["auto"] | list[tuple[str, str]] = _f(
        "auto", "homologous chromosome pairs: auto (from SYN alignments) or a list of [reference, query]"
    )
    min_syn_frac: float = _f(
        0.1,
        "keep a chromosome pair with at least this share of the best pair's SYN (1 = best only)",
        ge=0,
        le=1,
    )
    gapmax: int = _f(50_000, "maximum distance (bp) from a sequence to each anchor", ge=0)
    tol: int = _f(2_000, "minimum slack (bp) added to both ends of the projected window", ge=0)
    tol_frac: float = _f(0.1, "slack as a fraction of the distance to the farthest anchor", ge=0, le=1)
    flank: int = _f(1_000, "anchors closer than this (bp) are reported as flank_anchors evidence", ge=0)
    min_aligned_frac: float = _f(
        0.5, "an empty window counts as an aligned locus when at least this fraction is aligned", ge=0, le=1
    )

    @field_validator("types", mode="before")
    @classmethod
    def _split(cls, value: Any) -> Any:
        return value.split(",") if isinstance(value, str) else value

    @field_validator("types")
    @classmethod
    def _known(cls, value: list[str]) -> list[str]:
        types: list[str] = []
        for t in (v.strip().upper() for v in value):
            if t not in SYRI_CLASSES:
                raise ValueError(f"unknown SyRI class {t!r}, expected any of {', '.join(SYRI_CLASSES)}")
            if t not in types:
                types.append(t)
        if not types:
            raise ValueError("at least one anchor type is required")
        return types


class Validation(_Section):
    enabled: bool = _f(True, "confirm candidate pairs by aligning their flanks with BLASTN")
    mode: Literal["colinear", "coverage"] = _f(
        "colinear", "colinear: flanks must join the matching flank; coverage: any HSP counts (v3)"
    )
    reciprocal: bool = _f(True, "both loci of a pair must pass, not only the tested one")
    window: int = _f(10_000, "flank length (bp) taken on each side of both loci", ge=1)
    min_hsp: int = _f(200, "minimum HSP length (bp)", ge=1)
    min_frac: float = _f(0.2, "minimum aligned fraction of each flank", ge=0, le=1)
    min_bp: int = _f(1_000, "minimum aligned bases of each flank", ge=0)
    max_offset: int = _f(500, "how far (bp) a flank HSP may run past the other locus", ge=0)


class Dotplots(_Section):
    pairs: bool = _f(True, "a Gepard dotplot for every confirmed pair")
    genome: bool = _f(True, "a whole-genome dotplot of the SyRI alignments")


HEX_COLOR = r"^#[0-9A-Fa-f]{6}$"


def _color(default: str, text: str) -> Any:
    return Field(default, description=text, pattern=HEX_COLOR, json_schema_extra={"color": True})


class Figures(_Section):
    reference: str = _color("#3A3A3A", "reference species: chromosome lines and Venn circle")
    query: str = _color("#8C8C8C", "query species: chromosome lines and Venn circle")
    syntenic: str = _color("#DCDCDC", "syntenic blocks (SYN)")
    inversion: str = _color("#F4B97A", "inversions (INV)")
    translocation: str = _color("#B9B3E0", "translocations (TRANS, INVTR)")
    duplication: str = _color("#9CCBE8", "duplications (DUP, INVDP)")
    centromere: str = _color("#FCEFC4", "centromeres (from the BED files)")


class Species(_Section):
    label: str
    match: str
    genome: Path | None = None
    coords: Literal["header"] | Path = "header"
    header_regex: str = HEADER_REGEX
    lift: list[Path] = Field(default_factory=list)
    strip_suffix: str | None = None
    centromeres: Path | None = None
    min_identity: float = Field(0.99, ge=0.5, le=1)

    @field_validator("match")
    @classmethod
    def _match(cls, value: str) -> str:
        _compile(column_pattern(value))
        return value

    @property
    def column_regex(self) -> re.Pattern[str]:
        return re.compile(column_pattern(self.match))

    @field_validator("header_regex")
    @classmethod
    def _header(cls, value: str) -> str:
        missing = {"chr", "start", "end"} - set(_compile(value).groupindex)
        if missing:
            raise ValueError(f"missing named groups: {', '.join(sorted(missing))}")
        return value


class Run(_Section):
    label: str | None = None
    color: str | None = None
    unit: str | None = None
    path: Path | None = None


class OrthoFinderInput(_Section):
    root: Path
    runs: dict[str, Run] = Field(default_factory=dict)

    @field_validator("runs", mode="before")
    @classmethod
    def _names(cls, value: Any) -> Any:
        if isinstance(value, list):
            return {str(name): {} for name in value}
        if isinstance(value, dict):
            return {str(name): run or {} for name, run in value.items()}
        return value


class Config(_Section):
    orthofinder: OrthoFinderInput
    species: dict[str, Species]
    reference: str
    query: str
    outdir: Path = Path("spot_out")
    threads: int | Literal["auto"] = "auto"
    profile: str = "default"
    resume: bool = True
    syri_out: Path | None = None
    minimap2: Minimap2 = Field(default_factory=Minimap2)
    syri: Syri = Field(default_factory=Syri)
    blastn: Blastn = Field(default_factory=Blastn)
    gepard: Gepard = Field(default_factory=Gepard)
    anchors: Anchors = Field(default_factory=Anchors)
    validation: Validation = Field(default_factory=Validation)
    dotplots: Dotplots = Field(default_factory=Dotplots)
    figures: Figures = Field(default_factory=Figures)

    @field_validator("species")
    @classmethod
    def _two(cls, value: dict[str, Species]) -> dict[str, Species]:
        if len(value) != 2:
            raise ValueError(f"exactly two species are required, got {len(value)}")
        for key in value:
            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", key):
                raise ValueError(f"species key {key!r} must be alphanumeric")
        return value

    @field_validator("threads")
    @classmethod
    def _threads(cls, value: int | str) -> int | str:
        if isinstance(value, int) and value < 1:
            raise ValueError("threads must be at least 1, or auto")
        return value

    @property
    def n_threads(self) -> int:
        """Worker threads for this run (``auto`` = performance cores)."""
        from .resources import resolve

        return resolve(self.threads)

    @model_validator(mode="after")
    def _roles(self) -> Config:
        for role in ("reference", "query"):
            if getattr(self, role) not in self.species:
                raise ValueError(f"{role} {getattr(self, role)!r} is not one of {', '.join(self.species)}")
        if self.reference == self.query:
            raise ValueError("reference and query must be different species")
        return self


SECTIONS: dict[str, type[BaseModel]] = {
    "minimap2": Minimap2,
    "syri": Syri,
    "blastn": Blastn,
    "gepard": Gepard,
    "anchors": Anchors,
    "validation": Validation,
    "dotplots": Dotplots,
    "figures": Figures,
}


def column_pattern(match: str) -> str:
    """``match`` as a regular expression: a plain name such as ``MS`` means "columns starting with MS"."""
    return "^" + re.escape(match) if re.fullmatch(r"[A-Za-z0-9_.\-]+", match) else match


def _compile(pattern: str) -> re.Pattern[str]:
    try:
        return re.compile(pattern)
    except re.error as e:
        raise ValueError(f"invalid regex {pattern!r}: {e}") from None


def _read_yaml(path: Path) -> dict[str, Any]:
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except OSError as e:
        raise ConfigError(f"cannot read {path}: {e.strerror}") from None
    except yaml.YAMLError as e:
        raise ConfigError(f"{path}: {e}") from None
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: expected a mapping at the top level")
    return data


def merge(base: Mapping[str, Any], over: Mapping[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for key, value in over.items():
        if isinstance(value, Mapping) and isinstance(out.get(key), Mapping):
            out[key] = merge(out[key], value)
        else:
            out[key] = value
    return out


def assignment(text: str) -> dict[str, Any]:
    key, sep, value = text.partition("=")
    if not sep or not key.strip():
        raise ConfigError(f"expected KEY=VALUE, got {text!r}")
    node: Any = yaml.safe_load(value) if value.strip() else None
    for part in reversed(key.strip().split(".")):
        node = {part: node}
    out: dict[str, Any] = node
    return out


def builtin_profiles() -> dict[str, str]:
    folder = resources.files("spot").joinpath("profiles")
    return {
        entry.name.removesuffix(".yaml"): entry.read_text(encoding="utf-8")
        for entry in sorted(folder.iterdir(), key=lambda e: e.name)
        if entry.name.endswith(".yaml") and not entry.name.startswith(".")
    }


def load_profile(name: str, base: Path | None = None) -> tuple[str, dict[str, Any]]:
    if name.endswith((".yaml", ".yml")):
        path = Path(name).expanduser()
        if not path.is_absolute() and base is not None:
            path = base / path
        if not path.is_file():
            raise ConfigError(f"profile file not found: {path}")
        data = _read_yaml(path)
    else:
        profiles = builtin_profiles()
        if name not in profiles:
            raise ConfigError(f"unknown profile {name!r}; available: {', '.join(profiles)}")
        data = yaml.safe_load(profiles[name]) or {}
    description = str(data.pop("description", "")).strip()
    extra = sorted(set(data) - set(SECTIONS))
    if extra:
        allowed = ", ".join(SECTIONS)
        raise ConfigError(f"profile {name!r} sets {', '.join(extra)}; profiles only set {allowed}")
    return description, data


def profile_settings(name: str) -> dict[str, Any]:
    _, data = load_profile(name)
    return {
        section: model.model_validate(data.get(section, {})).model_dump(mode="json")
        for section, model in SECTIONS.items()
    }


def load(path: Path, *, profile: str | None = None, overrides: Mapping[str, Any] | None = None) -> Config:
    path = path.expanduser().resolve()
    user = _read_yaml(path)
    name = profile or user.get("profile") or "default"
    _, settings = load_profile(str(name), base=path.parent)
    data = merge(settings, user)
    if overrides:
        data = merge(data, overrides)
    data["profile"] = str(name)
    try:
        config = Config.model_validate(data)
    except ValidationError as e:
        raise ConfigError(_describe(e, path)) from None
    resolved: Config = _resolve(config, path.parent)
    return resolved


def dump(config: Config) -> str:
    return yaml.safe_dump(config.model_dump(mode="json"), sort_keys=False, allow_unicode=True)


def _describe(error: ValidationError, path: Path) -> str:
    lines = [f"invalid configuration in {path}:"]
    for item in error.errors():
        where = ".".join(str(p) for p in item["loc"]) or "(top level)"
        message = item["msg"].removeprefix("Value error, ")
        lines.append(f"  {where}: {message}")
    return "\n".join(lines)


def _resolve(value: Any, base: Path) -> Any:
    if isinstance(value, Path):
        p = value.expanduser()
        return p if p.is_absolute() else Path(os.path.normpath(base / p))
    if isinstance(value, BaseModel):
        fields = type(value).model_fields
        return value.model_copy(update={name: _resolve(getattr(value, name), base) for name in fields})
    if isinstance(value, list):
        return [_resolve(v, base) for v in value]
    if isinstance(value, dict):
        return {k: _resolve(v, base) for k, v in value.items()}
    return value
