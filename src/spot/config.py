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


class Minimap2(_ToolOptions):
    x: Literal["asm5", "asm10", "asm20"] = "asm10"


class Syri(_ToolOptions):
    nosnp: bool = True


class Blastn(_ToolOptions):
    evalue: float = 1e-10
    perc_identity: float = Field(80.0, ge=0, le=100)


class Gepard(_ToolOptions):
    matrix: str = "edna.mat"
    word: int = Field(10, ge=1)


class Anchors(_Section):
    types: list[str] = Field(default_factory=lambda: ["SYN", "INV"])
    same_chromosome: bool = True
    homology: Literal["auto"] | list[tuple[str, str]] = "auto"
    min_syn_frac: float = Field(0.1, ge=0, le=1)
    gapmax: int = Field(50_000, ge=0)
    tol: int = Field(2_000, ge=0)
    tol_frac: float = Field(0.1, ge=0, le=1)
    flank: int = Field(1_000, ge=0)
    min_aligned_frac: float = Field(0.5, ge=0, le=1)

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
    enabled: bool = True
    mode: Literal["colinear", "coverage"] = "colinear"
    reciprocal: bool = True
    window: int = Field(10_000, ge=1)
    min_hsp: int = Field(200, ge=1)
    min_frac: float = Field(0.2, ge=0, le=1)
    min_bp: int = Field(1_000, ge=0)
    max_offset: int = Field(500, ge=0)


class Dotplots(_Section):
    pairs: bool = True
    genome: bool = True


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
        _compile(value)
        return value

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
    threads: int = Field(4, ge=1)
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

    @field_validator("species")
    @classmethod
    def _two(cls, value: dict[str, Species]) -> dict[str, Species]:
        if len(value) != 2:
            raise ValueError(f"exactly two species are required, got {len(value)}")
        for key in value:
            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", key):
                raise ValueError(f"species key {key!r} must be alphanumeric")
        return value

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
}


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
