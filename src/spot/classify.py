from __future__ import annotations

from dataclasses import dataclass, field

from .orthofinder import OrthoFinderRun
from .positional import EVALUABLE, SYNTENIC_ORTHOLOG, Test

SYNTENIC = "syntenic_ortholog"
NOT_SYNTENIC = "evaluable_not_syntenic"
NOT_EVALUABLE = "not_evaluable"
UNASSIGNED = "unassigned"


@dataclass
class RunResult:
    run: OrthoFinderRun
    reference: str
    query: str
    tests: dict[str, Test]
    status: dict[str, str] = field(default_factory=dict)
    direction: dict[str, str] = field(default_factory=dict)
    pairs: dict[str, list[tuple[str, str]]] = field(default_factory=dict)

    def zone(self, gene: str) -> str:
        g = self.run.genes[gene]
        if not g.assigned:
            return f"{g.species}_unassigned"
        return self.og_zone(g.orthogroup)

    def og_zone(self, og: str) -> str:
        members = self.run.orthogroups[og].members
        present = [sp for sp in (self.reference, self.query) if members[sp]]
        return "shared" if len(present) == 2 else f"{present[0]}_only"

    def folder(self, og: str) -> str:
        zone = self.og_zone(og)
        return self.status[og] if zone == "shared" else zone


def classify(run: OrthoFinderRun, tests: dict[str, Test], reference: str, query: str) -> RunResult:
    result = RunResult(run, reference, query, tests)
    for og in run.orthogroups.values():
        if result.og_zone(og.id) != "shared":
            continue
        member_tests = [tests[g] for sp in (reference, query) for g in og.members[sp]]
        directions = sorted(
            {t.species for t in member_tests if t.klass == SYNTENIC_ORTHOLOG}, key=lambda sp: sp != reference
        )
        other = {reference: query, query: reference}
        result.direction[og.id] = "+".join(f"{sp}->{other[sp]}" for sp in directions)
        if directions:
            result.status[og.id] = SYNTENIC
        elif any(t.klass in EVALUABLE for t in member_tests):
            result.status[og.id] = NOT_SYNTENIC
        else:
            result.status[og.id] = NOT_EVALUABLE
        pairs = set()
        for t in member_tests:
            if t.klass == SYNTENIC_ORTHOLOG:
                for o in t.orthologs:
                    pairs.add((t.gene, o) if t.species == reference else (o, t.gene))
        result.pairs[og.id] = sorted(pairs)
    return result
