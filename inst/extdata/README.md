# Example data

`example.mrd` holds three precursor blocks copied verbatim out of a full
`miRBase.mrd`, chosen so that the example exercises every branch of the
classification:

| Precursor | Reads in mrd | Classified | Why |
|---|---|---|---|
| `hsa-mir-134` | 29 | 28 | covers all eight isomiR categories; the precursor used in `docs/algorithm.md` |
| `hsa-mir-320c-1` | 15 | 13 | a single annotated mature product (`M` arm), so the `-M` suffix rule is exercised |
| `hsa-mir-548aj-2` | 5 | 2 | most reads match no annotated arm and are dropped |

Running `bin/quantify_isomirs.R` on it yields 43 classified reads, stored in
`example.expected.tsv` and checked by `tests/test_example.R`.

Source: a small RNA-seq library of the SenomiR cell culture dataset, sequenced
in-house and quantified with miRDeep2 `quantifier.pl` against miRBase. The three
blocks are unmodified, and their classification is identical to the one obtained
when the full mrd file is processed.
