# isomiRScope

![isomiRScope](docs/isomiRScope.png)

isomiRScope classifies isomiRs from the `miRBase.mrd` file produced by the
miRDeep2 `quantifier.pl` script. For every collapsed read that miRDeep2 assigned
to a miRNA precursor it reports the arm the read belongs to, whether the read is
canonical or a 5' / 3' length variant, and whether it carries a mismatch against
the precursor.

The tool takes one mrd file and writes one table. It does not map reads, does
not need a genome and does not depend on the rest of the analysis pipeline it
was extracted from.

## Classification

| `isomir_name` | Meaning |
|---|---|
| `canonical` | both ends match the annotated mature sequence |
| `iso-5` | 5' end shifted |
| `iso-3` | 3' end shifted |
| `iso-3-5` | both ends shifted |
| `poly_*` | as above, plus at least one mismatch against the precursor |

The full rules, including the tolerance windows used to assign a read to an arm,
are in [docs/algorithm.md](docs/algorithm.md).

## Installation

```bash
conda env create -f environment.yml
conda activate isomirscope
```

## Usage

Produce the input with miRDeep2, for example:

```bash
quantifier.pl -p hairpin-hsa.fa -m mature-hsa.fa -r reads_collapsed.fa -y default
# writes expression_analyses/expression_analyses_default/miRBase.mrd
```

Then classify:

```bash
Rscript bin/quantify_isomirs.R \
  --readsFile=expression_analyses/expression_analyses_default/miRBase.mrd \
  --outFile=isomirs.tsv \
  --cores=4
```

| Argument | |
|---|---|
| `--readsFile` | `miRBase.mrd` written by `quantifier.pl` |
| `--outFile` | output path for the tab-separated isomiR table |
| `--cores` | cores used for classification, default 4 |

## Output

One row per collapsed read:

| Column | Meaning |
|---|---|
| `read_id` | collapsed read identifier from the mrd file |
| `count` | number of reads collapsed into that sequence |
| `value` | alignment string against the precursor |
| `seq_start`, `seq_end` | read boundaries in precursor coordinates |
| `mir_name` | precursor name with the assigned arm appended |
| `isomir_name` | `isomir_type` combined with `polymorphic` |
| `isomir_type` | `canonical`, `iso-5`, `iso-3` or `iso-3-5` |
| `polymorphic` | `poly` if the read carries a mismatch, otherwise `NA` |

Reads that cannot be assigned to an annotated arm are not reported.

## Example

`inst/extdata/example.mrd` is a small mrd file covering every isomiR type.
Running the tool on it reproduces `inst/extdata/example.expected.tsv`:

```bash
Rscript tests/test_example.R
```

## Citation

TODO: manuscript reference and Zenodo DOI once published.

## License

GPL-3. See [LICENSE](LICENSE) and [COPYRIGHT](COPYRIGHT).

isomiRScope was developed at [TAmiRNA GmbH](https://www.tamirna.com) as part of
the miND small RNA sequencing pipeline.
