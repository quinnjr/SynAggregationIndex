# SynAggregationIndex

Compute an α-synuclein aggregation index, `oligomer / (oligomer + monomer)`,
from separate per-sample quantifications of oligomeric and monomeric α-syn
(e.g. from oligomer-specific ELISAs, proximity-ligation assays, or
single-molecule array readouts).

Elevated oligomer:monomer ratios track Lewy-body pathology in PD, DLB, and
MSA research cohorts.

**Research tool, not a clinical diagnostic.**

## Input

Two sample-value CSVs (same format emitted by `SynProteinFilter`):

- `oligomer` — oligomeric α-syn quant per sample
- `monomer` — monomeric α-syn quant per sample

`log_inputs` defaults to `false`, i.e. inputs are assumed to be **linear/raw
values** — the format produced directly by assay CSVs such as the
Parkinsons pipeline's oligomer/monomer quantification files. This is the
default used by the real Parkinsons pipeline, which feeds raw linear assay
CSVs into this plugin.

If you instead chain this plugin directly from `SynProteinFilter`'s
**default** output (which is log2-normalized), you **must** set
`log_inputs<TAB>true` in the parameters file so the plugin reverses the
log2 transform before combining oligomer and monomer values. Leaving
`log_inputs` at its default `false` while feeding log2 values will silently
compute the aggregation index on log2 quantities instead of linear ones.

Example `parameters.synagg.txt` for chaining from `SynProteinFilter`'s
default log2 output:

```
oligomer	oligomer.synproteinfilter.csv
monomer	monomer.synproteinfilter.csv
log_inputs	true
zero_guard	1e-6
```

## Parameters (`parameters.synagg.txt`, tab-delimited)

| Key | Required | Default | Meaning |
|---|---|---|---|
| `oligomer` | yes | — | Path to oligomer quant CSV |
| `monomer` | yes | — | Path to monomer quant CSV |
| `log_inputs` | no | `false` | Assumes linear/raw inputs by default; set `true` to reverse log2 first (required when chaining from `SynProteinFilter`'s default log2 output) |
| `zero_guard` | no | `1e-6` | Minimum denominator |

## Output (TSV)

```
sample    agg_index    oligomer    monomer
PD_001    0.47         220         250
CTRL_001  0.08         14          170
```
