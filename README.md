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

## Parameters (`parameters.synagg.txt`, tab-delimited)

| Key | Required | Default | Meaning |
|---|---|---|---|
| `oligomer` | yes | — | Path to oligomer quant CSV |
| `monomer` | yes | — | Path to monomer quant CSV |
| `log_inputs` | no | `false` | Reverse log2 before combining |
| `zero_guard` | no | `1e-6` | Minimum denominator |

## Output (TSV)

```
sample    agg_index    oligomer    monomer
PD_001    0.47         220         250
CTRL_001  0.08         14          170
```
