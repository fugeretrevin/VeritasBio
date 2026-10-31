# Structure assignments per landscape

**Resolved 2026-10-31.** Every row was measured by aligning the candidate chain to the
landscape's wild-type FASTA and checking whether each library position is present in
coordinates. Resolution, identity and coverage are as measured, not as advertised.

This file exists because the landscape registry carries **no structure field**. Before
it, `spatial_filter.DEFAULT_PDB_PATH` silently fell back to avGFP (1EMA) whenever
`pdb_path` was `None`, applying one protein's contact shells, C-alpha distances and
burial estimates to a different target — with entirely plausible-looking output. That
fallback now warns (`_resolve_pdb_path`), but a warning only helps if there is a correct
answer to supply, which is what this table is.

## Coverage vocabulary

- **covered** — the position is in an exact-identity alignment block and the residue is
  observed in coordinates.
- **unmapped** — the position is not in an identical alignment block: either absent from
  the construct, or the construct carries a different residue there.
- **unobserved** — the residue is in the sequence but has no coordinates (a missing
  residue, commonly a terminus or a disordered loop).

## Assignments

| landscape | level | PDB | method / resolution | identity | library coverage | caveats |
| --- | --- | --- | --- | --- | --- | --- |
| `trpb` | development | **8EH0** | X-ray | — | in use | local copy; the only structure this project has scored against at length |
| `gb1` | development | **2QMT** | X-ray 1.05 Å | 1.00 | **4/4** | local copy present. `2GI9` (1.14 Å) equivalent. `7QTS`/`7QTR` are numbered +1 — do not use without renumbering |
| `her2` | development | **1N8Z** | X-ray | — | — | local copy; coverage against the 10-site library not re-measured |
| `cr9114_h1` | final_test | **5WL2** | X-ray 2.0 Å | 0.98 | **14/16** | 74 and 106 unmapped. `4FQH`/`4NZT`/`5CJQ`/`5CJS`/`4FQY`/`8UZP` are ~0.88 identity and 4/16 — **not** the corpus VH at the CDRs |
| `cr6261_h1`, `cr6261_h9` | final_test | **5WKZ** | X-ray 1.85 Å | 0.975 | **10/11** | 74 unmapped. Other candidates ~0.87 identity and 1/11 |
| `creilov` | final_test | **8QI9** | serial synchrotron, RT, 1.87 Å | 0.992 | **14/15** | CrPhotLOV1 dark state. 113 unobserved. Light-state sister `8QIB` (2.45 Å, 14/15); `8QIA` (2.50 Å, 13/15). `1N9L` is 0.908 identity. `2V0U`/`2V0W`/`2V1A` unusable |
| `eqfp611_blue`, `eqfp611_red`, `eqfp611_combined` | final_test | **3M24** | X-ray 2.20 Å | 0.983 | **10/13** | See the warning below. 63 and 174 **differ** from the corpus WT; 231 unobserved. Other family hits (`4NWS`, `4KPI`, `4KGE`, `8SXC`, `3BX9`) cover 1–3 of 13 |
| `hiv1_protease` | reserve | **1WBK** | X-ray 2.0 Å | 0.93–0.98 | **11/11** | `1WBM`, `2BQV`, `5T84` (1.65 Å) equivalent. Variant constructs |
| `rhla_rha`, `rhla_rls`, `rhla_pct_rha` | reserve | **8IK2** | X-ray 2.15 Å | 1.00 | **9/11** | 73 and 74 **unobserved**. Library positions taken from the `mutant` column, not diffed against WT |
| `sacas9` | reserve | **5AXW** | X-ray 2.7 Å | 0.995 | **8/8** | `5CZZ` (2.6 Å), `7EL1` (2.22 Å) equivalent. `8ZCY`/`8ZCZ`/`8ZD0`/`8ZDA` are 5/8 — 985, 986, 991 unmapped |
| `spcas9` | reserve | **5F9R** | X-ray 3.4 Å | ≥0.999 | **8/8** | `5VZL`, `6MCC`, `8KAH` equivalent. Each has 3–5 polymer entities (guide RNA / DNA): none is apo |

## Two warnings that change how a structure may be used

**eqFP611 carries the wrong residue at two library positions.** The corpus WT is
`eqFP611(mTagBFP2).fasta`, and mTagBFP2 is heavily engineered from eqFP611. 3M24 is
"TagBFP" and matches at 0.983 overall, but **63 and 174 differ from the corpus WT**. A
contact shell or burial estimate computed at those two sites describes a different
residue than the one the library mutates. Either exclude 63 and 174 from any
structure-derived term, or model the substitution explicitly. Do not treat 10/13 as if
it were 13/13 minus three missing atoms — three of the thirteen are *wrong*, not absent.
This is the same failure the avGFP fallback produced, in miniature.

**Apo versus complex is unresolved for every entry.** Non-polymer entity *counts* were
measured; the ligands were not named. 3M24 has four non-polymer entities and 8QI9 has
two, and in both cases a chromophore is the obvious candidate — LOV domains carry FMN —
but neither is confirmed. All SpCas9 entries are complexes with nucleic acid. Any term
that depends on the binding site being empty, or on a cofactor being present, is
currently unjustified.

## Not verified

- Ligand identity for any entry, hence apo/complex status.
- HER2 coverage against its 10-site library.
- Whether the three eqFP611 endpoints and the three RhlA endpoints should share one
  structure assignment — they share a source table, so presumably yes, but the library
  positions were derived once per source file, not once per endpoint.
- Resolution reported as 0.00 for `1Q10`, `1GB4`, `1MPE`, `7QTS`, `7QTR`, most likely
  NMR or no recorded resolution. Not checked against the entry's method field.

## Provenance

Structures were checked with `gemmi` against `files.rcsb.org/download/<ID>.cif`; chains
matched by `label_seq`, then aligned to the wild-type FASTA from
`downloads/combingym/Data/WT Sequence/`. Library positions were derived by diffing each
clean table's genotypes against that FASTA, except RhlA which used the `mutant` column.
The corpus is pinned at CombinGym `c5dccd514a9acfe1089a5f1555faea98afd78019` and
verifies 29/29.
