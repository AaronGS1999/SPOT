# Output

`MS` and `MO` stand for the keys of the reference and query species; `NUMT_2` for an OrthoFinder run.

```
outdir/
├── spot.log                  what was done, with every command run
├── run.yaml                  settings in force
├── versions.tsv              version and path of every program
├── homology.tsv              homologous chromosome pairs
├── calibration.tsv           hits observed and expected by chance, per run
├── venn_panel.pdf|png|svg    Venn diagrams of all runs side by side
├── synteny.pdf|png           homologous chromosomes with the confirmed pairs
├── genome/
│   ├── dotplot.pdf|png       all chromosomes against all
│   └── genomes.paf, syri.out, …   only when SPOT aligned the genomes
└── NUMT_2/                   one folder per OrthoFinder run
    ├── summary.txt
    ├── orthogroups.tsv
    ├── pairs.tsv
    ├── sequences.tsv
    ├── venn.pdf|png|svg
    ├── syntenic_ortholog/OG…/
    ├── evaluable_not_syntenic/OG…/
    ├── not_evaluable/OG…/
    ├── MS_only/OG…/, MO_only/OG…/
    └── unassigned/MS.fa, MO.fa
```

Each orthogroup folder holds:

- `OG….fa`: all its sequences;
- `OG….pairs.fa`: the sequences in confirmed pairs (syntenic orthogroups only);
- `OG…_tree.txt`: its OrthoFinder gene tree (Newick), when OrthoFinder built one (4 or more sequences);
- `dotplots/<MS id>__<MO id>.png`: a dotplot of each confirmed pair (see figures).

All coordinates are 1-based and inclusive. Tables are tab-separated with a header line.

## Tables

### `orthogroups.tsv`: one row per orthogroup

| Column | Content |
|---|---|
| `orthogroup` | OrthoFinder orthogroup ID |
| `zone` | `shared`, or `MS_only` / `MO_only` when only one species has members |
| `status` | shared orthogroups only: `syntenic_ortholog`, `evaluable_not_syntenic` or `not_evaluable` ([method](method.md#7-classify-the-orthogroups)) |
| `direction` | species from which the confirmed pairs were found: `MS->MO`, `MO->MS`, or both joined by `+` |
| `n_MS`, `n_MO` | number of sequences of each species |
| `tree` | `yes` if OrthoFinder built a gene tree |
| `q_value` | smallest q-value of its members' positional hits |
| `tree_orthologs` | confirmed pairs also called orthologous by the OrthoFinder gene trees, as `k/n` |
| `pairs` | confirmed pairs, `MS id|MO id`, separated by `;` |

### `pairs.tsv`: one row per confirmed pair

| Column | Content |
|---|---|
| `orthogroup` | orthogroup of the pair |
| `MS_id`, `MS_chr`, `MS_start`, `MS_end`, `MS_strand` | reference sequence and its position |
| `MO_id`, `MO_chr`, `MO_start`, `MO_end`, `MO_strand` | query sequence and its position (after lifting) |
| `found_from` | test that found the pair: `MS->MO`, `MO->MS`, or both joined by `+` |
| `p_value`, `q_value` | of the positional hit (the smaller one when found from both sides) |
| `tree_orthology` | `ortholog` if the OrthoFinder gene trees also call the pair orthologous, `not_ortholog` if not, empty when OrthoFinder wrote no orthologue table |

Strand is the orientation of the sequence, as given to OrthoFinder, relative to the genome.

### `sequences.tsv`: one row per sequence

| Column | Content |
|---|---|
| `id` | sequence ID |
| `species` | species key |
| `orthogroup` | orthogroup, empty if unassigned |
| `zone` | `shared`, `MS_only`, `MO_only`, `MS_unassigned` or `MO_unassigned` |
| `og_status`, `og_direction` | `status` and `direction` of its orthogroup (shared only) |
| `class` | result of its positional test ([classes](method.md#4-test-the-position-of-every-sequence)) |
| `evidence` | `flank_anchors` when both anchors are within `flank` (1 kb) of the sequence, `bracketing_anchors` when further away (up to `gapmax`) |
| `block_type` | `collinear` when both anchors are SYN or INV, `rearranged` otherwise; `_two_blocks` when they belong to two different SyRI blocks |
| `syri_blocks` | SyRI block IDs of the anchors |
| `chr`, `start`, `end`, `strand` | position in its genome |
| `location` | `located`, or why it was not: `no_coordinates`, `lift_failed`, `off_assembly`, `out_of_range`, `sequence_mismatch` |
| `identity` | identity between the sequence and the genome at its position (1 = identical) |
| `anchor_gap_left`, `anchor_gap_right` | distance (bp) to the nearest anchor on each side |
| `other_chr`, `window_start`, `window_end` | projected window on the other genome |
| `window_aligned_frac` | for empty windows, fraction covered by SyRI alignments |
| `orthologs_in_window` | members of the same orthogroup in the window, confirmed by their flanks |
| `others_in_window` | other sequences of the run in the window |
| `unconfirmed` | members of the same orthogroup in the window that were not confirmed: their flanks did not align, or the hit was not significant |
| `p_value` | probability of finding a member of the orthogroup in a window this long by chance ([method](method.md#5-test-the-significance-of-each-hit)) |
| `q_value` | Benjamini-Hochberg q-value of the hit within its run |
| `flank_left_bp`, `flank_right_bp` | colinear aligned bases of this sequence's left and right flank (best candidate) |
| `flank_back_left_bp`, `flank_back_right_bp` | the same for the candidate's flanks (reciprocal check) |
| `flank_orientation` | `direct`, or `inverted` when the pair lies inside an inversion |
| `description` | rest of the FASTA header |

### `homology.tsv`: one row per homologous chromosome pair

| Column | Content |
|---|---|
| `reference_chr`, `query_chr` | the pair |
| `syn_bp` | bases of the reference chromosome in syntenic alignments with this query chromosome |
| `fraction_of_reference_syn` | share of all syntenic bases of the reference chromosome |

### `calibration.tsv`: one row per run

| Column | Content |
|---|---|
| `run` | OrthoFinder run |
| `observed_hits` | windows holding a member of the same orthogroup, before significance and flanks |
| `permutations` | orthogroup-label permutations |
| `null_mean`, `null_sd`, `null_p95` | hits after shuffling the labels: mean, standard deviation and 95th percentile |
| `expected_false` | `null_mean / observed_hits`: share of the hits expected by chance |

### `summary.txt`

Counts per run: orthogroups and sequences per zone, shared orthogroups by status and direction, confirmed
pairs, every sequence by zone and class, positional hits and how many are significant, agreement of the
confirmed pairs with the gene trees, the calibration, and the syntenic tests by evidence and block type.

### `versions.tsv`

Program, version and path of SPOT, Python, minimap2, SyRI, BLASTN, Java, Gepard and the Python libraries
(including SciPy and Biopython).

## Figures

### Venn diagram (`venn`, `venn_panel`)

One diagram per run. Each circle is a species, labelled with the number of its sequences in the run. The
outer part of each circle counts the orthogroups with members of that species only, their sequences and,
below, the unassigned sequences. The overlap holds the shared orthogroups, split in two: syntenic
orthologs on top, in the colour of the run, and the other shared orthogroups below, in a pale shade
(evaluable but not syntenic, or not evaluable; `orthogroups.tsv` tells them apart). Each part gives the
number of orthogroups and of sequences of each species. `venn_panel` places the runs side by side.

### Synteny (`synteny`)

One panel per pair of homologous chromosomes, reference above (dark line) and query below (grey line), to
scale in Mbp. Ribbons join the SyRI blocks: syntenic (grey), inversions (orange), translocations (purple)
and duplications (blue). Pale yellow boxes are centromeres, when given. Curves join the confirmed pairs,
one colour per run, and dots on the chromosomes mark the sequences classed `syntenic_ortholog`. The small
area plots above the reference and below the query show, for each run, the number of sequences per 250 kb
along the chromosome, scaled to the run's maximum.

### Whole-genome dotplot (`genome/dotplot`)

All reference chromosomes (x axis) against all query chromosomes (y axis), each axis to scale. Each SyRI
alignment is a short segment coloured by type, so syntenic regions form diagonals, inversions run against
them and translocations and duplications sit off the diagonal. Dots are confirmed pairs, one colour per run.

### Pair dotplots (`dotplots/*.png`)

Gepard dotplot (Krumsiek et al., 2007) of the two loci of a confirmed pair with their flanks (10 kb each
side). A continuous diagonal through the centre means the flanks are colinear; a gap in the diagonal at
the centre is the insertion itself or a local difference.

All colours can be changed in **Parameters › Figure colours** (species, block types, centromeres) and in the
run table of **Project** (one colour per run).
