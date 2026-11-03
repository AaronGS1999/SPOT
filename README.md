# SPOT

Synteny-based Positional Orthology Test. Takes OrthoFinder orthogroups built from sequences of any kind
(not only genes) and checks, with SyRI on a whole-genome alignment of two species, whether the members of
each orthogroup sit at orthologous positions of homologous chromosomes.

## Install

macOS (Apple Silicon or Intel) and Linux.

```sh
mamba env create -f environment.yml
mamba activate spot
pip install -e .
spot check
```

On Apple Silicon the environment solves natively. If the solver complains about the macOS version, prefix the
first command with `CONDA_OVERRIDE_OSX=14.0`.

## Desktop app

Build the app once with `python packaging/build.py macos` (gives `dist/SPOT.app`, for Apple Silicon and
Intel) or `python packaging/build.py linux` (gives `dist/SPOT-linux/`, run its `install.sh`). Open it like any
other app: the first launch sets up its own environment with Python, Qt, BLAST, minimap2, SyRI, Gepard and
Java (a few minutes, needs internet; kept in `~/Library/Application Support/SPOT` or
`~/.local/share/spot`), later launches start at once. From a development install, `spot gui` opens the same
app.

The app edits the same YAML configuration as the command line. **Project** sets the OrthoFinder folder
(runs are found automatically; choose which to use, their labels, colours and units), the two species and
the output folder. **Parameters** shows every engine setting with its value and whether it comes from the
profile or from the project; only project values are saved. **Run** checks the configuration, runs the
engine in its own process (with a progress bar per stage and Cancel) and keeps a copy of exactly what ran
in `outdir/spot_app.yaml`. **Results** shows the summary, filterable tables, the figures and, for each
confirmed pair, its Gepard dotplot and gene tree.

## Usage

```sh
spot run     -c config.yaml        # whole analysis
spot config  -c config.yaml        # resolved configuration
spot inspect -c config.yaml        # what SPOT reads from each OrthoFinder run
spot check   -c config.yaml        # external tools needed by this configuration
spot profiles                      # built-in profiles
spot profiles v3
```

See `examples/config.yaml`. Settings are applied in this order, each overriding the previous one:
built-in defaults, profile, configuration file, command line.

```sh
spot run -c config.yaml -p relaxed -x asm20 --anchors SYN,INV,TRANS --set validation.min_frac=0.3
```

Options of external tools use the tools' own names (`minimap2: {x: asm10}`, `blastn: {evalue: 1e-10}`,
`syri: {invgaplen: 1000}`). Anchor types are SyRI classes: `SYN`, `INV`, `TRANS`, `INVTR`, `DUP`, `INVDP`.

| profile | anchors | flank validation |
|---|---|---|
| `default` | SYN, INV | aligned fraction and minimum length |
| `relaxed` | all six, homologous chromosome only | as default |
| `v3` | all six, homologous chromosome only | >= 1 kb per flank |

### Input

- `orthofinder.root`: a folder with one or more OrthoFinder results. Each first-level subfolder is a run;
  the latest `Results_*` inside it is used. `runs` restricts and labels them.
- Sequence coordinates come from the FASTA headers (`chr:start-end`, see `header_regex`) or from a table
  given in `coords` (columns `id chr start end`). `lift` converts them through AGP files or split tables
  when they refer to another assembly. Every sequence is checked against the genome before it is used.
- `syri_out`: an existing `syri.out` for the same two genomes. Without it SPOT runs
  `minimap2 -c --eqx` and `syri` and keeps the result in `outdir/genome` for later runs.

### Method

1. **Location.** Each sequence is placed on its genome from the FASTA header or a coordinate table, lifted
   through AGP/split tables when needed, and checked against the genome on both strands (exact, or at least
   `min_identity`, default 0.99). Sequences on scaffolds missing from the genome are `off_assembly`.
2. **Homology.** Chromosome pairs joined by SyRI SYN alignments; a pair is kept when it carries at least
   `anchors.min_syn_frac` (0.1) of the SYN of the best pair of either chromosome, so fusions and fissions
   give N:M pairs (`homology.tsv`). `anchors.homology` sets the pairs by hand.
3. **Positional test.** The nearest SyRI alignments of the chosen types (`anchors.types`, SYN and INV by
   default) on both sides of a sequence are projected onto the homologous chromosome. The window is widened
   by `max(tol, tol_frac * distance to the farthest anchor)` and searched for members of the same orthogroup.
   An empty window is `empty_locus_aligned` when at least `min_aligned_frac` of it is aligned by SyRI (the
   locus exists and holds no sequence) and `empty_locus_gap` otherwise (not evaluable).
4. **Flank validation.** Each candidate pair is aligned with BLASTN with `validation.window` bp of flank. In
   `colinear` mode a flank only counts where it joins the matching flank of the other locus (left to left on
   the plus strand, or left to right on the minus strand), and with `reciprocal` both loci must pass.
   `coverage` mode (profile `v3`) counts any HSP, as the v3 script did.

A shared orthogroup is a `syntenic_ortholog` when at least one pair is confirmed, `evaluable_not_syntenic`
when some member could be tested but none was confirmed, and `not_evaluable` otherwise.

### Threads

`threads: auto` (the default) uses the performance cores of the computer (on Apple Silicon the efficiency
cores are left out). BLASTN runs as one single-threaded search per worker over a shared database of all
flanked loci, each restricted to the candidates of its own loci; Gepard draws all dotplots of a worker in a
single Java process, which is about 15 times faster than one JVM per dotplot. Results are written with one
thread per OrthoFinder run, which helps on network drives.

Each stage keeps a checkpoint in `outdir/.spot/`, keyed by a hash of its inputs and settings, so a repeated
or interrupted run only redoes what changed (`resume: false` turns this off). Input digests are cached in
`~/.cache/spot` (`SPOT_CACHE`).

### Output

```
outdir/
  spot.log  run.yaml  versions.tsv  homology.tsv
  venn_panel.{pdf,png,svg}  synteny.{pdf,png}
  genome/dotplot.{pdf,png}            (and minimap2/SyRI files when SPOT aligned the genomes)
  <run>/
    orthogroups.tsv  pairs.tsv  sequences.tsv  summary.txt  venn.{pdf,png,svg}
    syntenic_ortholog/<OG>/           <OG>.fa  <OG>.pairs.fa  <OG>_tree.txt  dotplots/*.png
    evaluable_not_syntenic/<OG>/
    not_evaluable/<OG>/
    <SP>_only/<OG>/
    unassigned/<SP>.fa
```

`sequences.tsv` has one row per sequence with its location and identity to the genome, test class, anchors,
the aligned fraction of the projected window and the colinear flank bases in both directions;
`pairs.tsv` lists the confirmed pairs and the direction in which they were found.

## Development

```sh
pip install -e ".[test]"
pytest                 # tests (BLAST tests are skipped when blastn is not on PATH)
ruff check src tests   # lint
mypy src               # strict typing
```

Your own data never goes into the repository: keep it outside `code/` (for example in `../DATA`) and point the
configuration at it.
