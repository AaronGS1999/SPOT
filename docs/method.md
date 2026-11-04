# How SPOT works

SPOT answers one question for every orthogroup that two species share: do its members lie at the same
position in both genomes? A shared orthogroup only says that sequences are similar. Being in the same
place, between the same conserved neighbours, is what makes two copies positional orthologs, that is,
descendants of a single insertion or gene present in the common ancestor.

The analysis has six steps. Names in `code` are settings (see `examples/config.yaml`).

## 1. Read the orthogroups

For each OrthoFinder run SPOT reads `Orthogroups.tsv`, `Orthogroups_UnassignedGenes.tsv`, the sequences
and the gene trees. Each column of `Orthogroups.tsv` is assigned to one species by its name
(`species.<key>.match`, e.g. `MS` for `MS_NUMT_evento2`). Every orthogroup then falls in one zone:

- **shared**: it has members in both species;
- **only one species** (`MS_only`): all its members are from one species;
- **unassigned**: sequences OrthoFinder left out of every orthogroup.

Only shared orthogroups can be tested.

## 2. Place every sequence on its genome

The position of each sequence is read from its FASTA header (`chr1:2261774-2261923`, 1-based and inclusive)
or from a table. When the coordinates were taken on another assembly, they are converted through the lift
files in order: split tables (one chromosome cut into pieces) and AGP files (pieces placed into
chromosomes, for example by RagTag).

The genome sequence at those coordinates is then compared with the sequence itself, on both strands. A
sequence is **located** when they are identical or at least `min_identity` similar (0.99). Otherwise it is
left out with the reason: no coordinates, outside the lift files, on a scaffold missing from the genome
(`off_assembly`), outside its chromosome, or a sequence mismatch. This check guards against wrong or
shifted coordinates, the most common silent error in this kind of analysis.

## 3. Align the two genomes

SPOT needs a whole-genome alignment of the two species classified by SyRI (Goel et al., 2019): which
aligned stretches are syntenic (SYN), inverted (INV), translocated (TRANS, INVTR) or duplicated (DUP,
INVDP). Either an existing `syri.out` is given, or SPOT aligns the genomes with minimap2 (Li, 2018; preset
`asm10`, about 1 % divergence) and runs SyRI on the result.

From the syntenic alignments SPOT pairs the homologous chromosomes: a reference and a query chromosome are
homologous when the syntenic sequence joining them is at least `min_syn_frac` (10 %) of that of the best
partner of either chromosome. This keeps one-to-one pairs and also the many-to-one pairs left by
chromosome fusions or fissions (`homology.tsv`). The pairs can also be given by hand.

## 4. Test the position of every sequence

For each located sequence of a shared orthogroup, in both directions (reference to query and query to
reference):

1. **Anchors.** The nearest SyRI alignments on its left and on its right are taken as anchors, within
   `gapmax` (50 kb). Only the alignment types in `anchors.types` are used, SYN and INV by default, and only
   alignments joining homologous chromosomes. Translocations and duplications are left out by default
   because their anchors may point to a copy rather than to the orthologous locus.
2. **Projection.** Each anchor is projected onto the other genome by linear interpolation inside the
   alignment. The two projections bound the expected position of the orthologous locus.
3. **Window.** That interval is widened on both sides by `tol` (2 kb) or by `tol_frac` (10 %) of the
   distance to the farthest anchor, whichever is larger, to absorb indels between the anchors and the locus.
4. **Search.** The window is searched for sequences of the same OrthoFinder run.

Each sequence ends in one class:

| Class | Meaning | Evaluable |
|---|---|---|
| `syntenic_ortholog` | a member of the same orthogroup is in the window and the flanks confirm it (step 5) | yes |
| `syntenic_locus_flanks_not_confirmed` | a member of the same orthogroup is in the window, but the flanks do not align | yes |
| `syntenic_locus_other` | the window holds sequences of the run, none of the same orthogroup | yes |
| `empty_locus_aligned` | the window holds no sequence and at least `min_aligned_frac` (50 %) of it is aligned: the locus exists and is empty | yes |
| `empty_locus_gap` | the window holds no sequence but lies mostly in unaligned sequence: absence cannot be told | no |
| `no_anchor` | no anchor on one side or both | no |
| `anchors_too_far` | the nearest anchors are more than `gapmax` away | no |
| `anchors_split` | the two anchors point to different chromosomes | no |
| `anchors_inconsistent` | the projected anchors are much further apart than the anchors themselves (more than 3 times their span plus 10 kb) | no |
| `unplaced` | the sequence could not be located (step 2) | no |

## 5. Confirm each candidate pair by its flanks

Two sequences in the same window could still be independent insertions that happen to be close. Each
candidate pair is therefore checked with BLASTN (Camacho et al., 2009). Both loci are taken with `window`
(10 kb) of flank on each side and aligned. Only colinear alignments count: the left flank of one locus
must align to the left flank of the other and the right to the right on the same strand, or crosswise on
the opposite strand when the pair lies inside an inversion. An alignment that joins a flank to the wrong
side of the other locus, as repeats often do, is ignored. Alignments must have at least `perc_identity`
(80 %) identity and `min_hsp` (200 bp) length.

A pair is confirmed when both flanks have at least `min_bp` (1 kb) or `min_frac` (20 %) of their length
aligned, whichever is larger. With `reciprocal` (default) the flanks of both loci must pass. The
candidate becomes `syntenic_ortholog`; otherwise `syntenic_locus_flanks_not_confirmed`.

## 6. Classify the orthogroups

Each shared orthogroup takes the best result of its members:

- **syntenic_ortholog**: at least one confirmed pair, found from either species;
- **evaluable_not_syntenic**: at least one member could be tested, none confirmed;
- **not_evaluable**: no member could be tested.

The direction column says from which species the pairs were found (`MS->MO`, `MO->MS` or both). Pairs
found in both directions are the most robust.

## Reproducibility

- Every run writes the settings in force (`run.yaml`), the version of every program (`versions.tsv`) and a
  log (`spot.log`).
- The results depend only on the inputs and settings. There is no random step, the BLASTN and Gepard work
  is split in a fixed order, and each OrthoFinder run is tested only against its own sequences, so adding
  or removing a run does not change the others.
- Finished stages are kept in `outdir/.spot`, keyed by a SHA-256 of their inputs and settings. A repeated
  or interrupted run only redoes what changed.

## Performance

`threads: auto` uses the performance cores of the computer. BLASTN runs as single-threaded searches over a
shared database of all flanked loci, each limited to the candidates of its own loci, and Gepard draws all
the dotplots of a worker in a single Java process. On a laptop a run of three OrthoFinder sets (about
9,000 sequences) with an existing `syri.out` takes under a minute.
