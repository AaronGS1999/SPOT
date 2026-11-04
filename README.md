<p align="center">
  <img src="src/spot/gui/assets/spot.png" width="128" alt="SPOT">
</p>

<h1 align="center">SPOT</h1>

<p align="center">Synteny-based Positional Orthology Test</p>

SPOT tells which OrthoFinder orthogroups shared by two species sit at the same place in both genomes.
For every sequence it finds the syntenic anchors on each side, projects them onto the homologous
chromosome of the other species, looks for a member of the same orthogroup there, tests whether finding
it is more than chance (Benjamini-Hochberg FDR) and confirms the pair by aligning the flanks of both loci.
Each pair is also compared with the orthologues that OrthoFinder infers from its gene trees. It works with any sequence that has genomic coordinates: nuclear
insertions of organelle DNA (NUMTs, NUPTs), transposable elements, genes or non-coding RNAs.

- [How it works](docs/method.md), step by step
- [Output](docs/output.md): every table, column and figure
- [References](#references)

## Install

**App** (macOS on Apple Silicon or Intel, Linux on x86-64 or ARM)

```sh
python packaging/build.py macos     # dist/SPOT.app
python packaging/build.py linux     # dist/SPOT-linux/, then run its install.sh
```

Open SPOT like any other app. The first launch installs its tools (Python, Qt, BLAST+, minimap2, SyRI,
Gepard and Java) in `~/.spot` (macOS) or `~/.local/share/spot` (Linux); it takes a few minutes and needs an
internet connection. Later launches start at once.

**Command line**

```sh
conda env create -f environment.yml
conda activate spot
pip install -e .
spot check
```

## Use

1. **Project.** Choose the folder that holds the OrthoFinder runs and press *Find runs*. Set each
   species: its genome, and lift files if the coordinates in the sequence names were taken on another
   assembly. The panel under each species says how many sequences fall on its genome.
2. **Parameters.** The default profile is a sound start. Values changed in the project are marked and can
   be reset to the profile.
3. **Run.** *Check configuration*, then *Run SPOT*.
4. **Results.** Tables, figures, and the dotplot and gene tree of each confirmed pair.

Hover over any field for an explanation. The same project runs from the terminal:

```sh
spot run -c project.yaml                     # whole analysis
spot run -c project.yaml -p relaxed -t 8     # another profile, 8 threads
spot config -c project.yaml                  # settings in force
spot inspect -c project.yaml                 # what is read from each OrthoFinder run
spot profiles                                # built-in profiles
```

`examples/config.yaml` lists every setting. Settings apply in this order, each overriding the previous
one: built-in defaults, profile, configuration file, command line (`--set section.key=value`).

| Profile | Anchors |
|---|---|
| `default` | syntenic blocks and inversions (SYN, INV) |
| `relaxed` | also translocations and duplications within homologous chromosomes |

## Input

- **OrthoFinder results**, one run per subfolder (`orthofinder -d` for nucleotide sequences). The latest
  `Results_*` of each run is used.
- **Coordinates** of each sequence, read from its FASTA header (`chr1:120001-120150`, detected
  automatically) or from a table (`id chr start end`).
- **Two genome FASTA files**, the assemblies the coordinates refer to (or are lifted to).
- Optional: a `syri.out` of the two genomes, and centromere BED files for the synteny figure.

## Development

```sh
pip install -e ".[test,gui]"
pytest
ruff check src tests
mypy src
```

## Citation

If you use SPOT, please cite this repository and the tools it relies on (`CITATIONS.bib`).

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach
to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological), 57*(1), 289–300.
https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009).
BLAST+: Architecture and applications. *BMC Bioinformatics, 10*, 421.
https://doi.org/10.1186/1471-2105-10-421

Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T.,
Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: Freely available Python tools for
computational molecular biology and bioinformatics. *Bioinformatics, 25*(11), 1422–1423.
https://doi.org/10.1093/bioinformatics/btp163

Emms, D. M., & Kelly, S. (2019). OrthoFinder: Phylogenetic orthology inference for comparative genomics.
*Genome Biology, 20*, 238. https://doi.org/10.1186/s13059-019-1832-y

Goel, M., Sun, H., Jiao, W.-B., & Schneeberger, K. (2019). SyRI: Finding genomic rearrangements and local
sequence differences from whole-genome assemblies. *Genome Biology, 20*, 277.
https://doi.org/10.1186/s13059-019-1911-0

Grüning, B., Dale, R., Sjödin, A., Chapman, B. A., Rowe, J., Tomkins-Tinch, C. H., Valieris, R., Köster, J.,
& The Bioconda Team. (2018). Bioconda: Sustainable and comprehensive software distribution for the life
sciences. *Nature Methods, 15*(7), 475–476. https://doi.org/10.1038/s41592-018-0046-7

Harris, C. R., Millman, K. J., van der Walt, S. J., Gommers, R., Virtanen, P., Cournapeau, D., Wieser, E.,
Taylor, J., Berg, S., Smith, N. J., Kern, R., Picus, M., Hoyer, S., van Kerkwijk, M. H., Brett, M., Haldane,
A., del Río, J. F., Wiebe, M., Peterson, P., … Oliphant, T. E. (2020). Array programming with NumPy.
*Nature, 585*(7825), 357–362. https://doi.org/10.1038/s41586-020-2649-2

Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering, 9*(3),
90–95. https://doi.org/10.1109/MCSE.2007.55

Krumsiek, J., Arnold, R., & Rattei, T. (2007). Gepard: A rapid and sensitive tool for creating dotplots on
genome scale. *Bioinformatics, 23*(8), 1026–1028. https://doi.org/10.1093/bioinformatics/btm039

Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. *Bioinformatics, 34*(18), 3094–3100.
https://doi.org/10.1093/bioinformatics/bty191

Li, H. (2021). New strategies to improve minimap2 alignment accuracy. *Bioinformatics, 37*(23), 4572–4574.
https://doi.org/10.1093/bioinformatics/btab705

Pydantic Services Inc. (2024). *Pydantic* [Computer software]. https://docs.pydantic.dev

Shirley, M. D., Ma, Z., Pedersen, B. S., & Wheelan, S. J. (2015). *Efficient "pythonic" access to FASTA files
using pyfaidx* (PeerJ PrePrints 3:e970v1). https://doi.org/10.7287/peerj.preprints.970v1

The Qt Company. (2024). *Qt for Python (PySide6)* [Computer software]. https://doc.qt.io/qtforpython-6/

Van Rossum, G., & Drake, F. L. (2009). *Python 3 reference manual*. CreateSpace.

Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., Burovski, E., Peterson,
P., Weckesser, W., Bright, J., van der Walt, S. J., Brett, M., Wilson, J., Millman, K. J., Mayorov, N.,
Nelson, A. R. J., Jones, E., Kern, R., Larson, E., … SciPy 1.0 Contributors. (2020). SciPy 1.0: Fundamental
algorithms for scientific computing in Python. *Nature Methods, 17*(3), 261–272.
https://doi.org/10.1038/s41592-019-0686-2
