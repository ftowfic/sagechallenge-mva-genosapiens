# Track 2 methods

This document describes released research software.
It does not claim a submitted medication nomination or treatment benefit.
`RELEASE.json` identifies the selected source and exported files.
The [shared methods](portable-methods.md) identify disclosure, performance, and review limits.

## Abstract

Track 2 connects retained genotype evidence to explicit mechanism hypotheses and medication research nominations.
Inputs include diagnosis-guided variants or reviewed published genotypes.
The workflow preserves allele, phase, source, and mechanism uncertainty.
Declared readiness rules separate accepted mechanisms from unresolved or deferred hypotheses.
Public discovery uses pinned target, drug, regulatory, and LINCS resources.

Native tissue analysis compares genetic perturbation and drug signatures within eligible cell and tissue contexts.
It preserves overexpression, knockdown, and paired-method distinctions.
Reports retain conflicting conditions, missing evidence, failed work, and applicability limits.
Declared strength calculations and frozen-background calibration support relative research prioritization.
The report combines mechanism rationale, medication evidence, figures, complete supplements, and source references.
Finalization checks artifacts while preserving scientific incompleteness and actual review status.

Separate reporting software can regenerate supported retained results without original case inputs or repeated discovery.
Producer admission and direct scientific dependency checks constrain historical transitions.
Unsupported transitions preserve original packages, checkpoints, reports, and review records.

Explicit evidence provenance and separate completion states support review and reproducibility.
Cell-line relevance, incomplete target coverage, mechanistic uncertainty, and absent clinical validation limit interpretation.
Research nominations do not establish efficacy, suitable exposure, safety, or treatment eligibility.

## Inputs and mechanism characterization

Diagnosis-guided execution accepts a VCF with a diagnosis identifier and optional HPO context.
Public reconstruction accepts an authenticated published-genotype bundle and actual input and plan reviews.
Shared retention rules preserve supported alleles before genotype and mechanism construction.
An uncertain phase or unresolved partner remains explicit.

Mechanism assessment separates loss of function, gain of function, and unsupported or uncertain effects.
It preserves the asserted molecular effect, pathway rationale, downstream consequence, and supporting evidence.
Mechanism readiness does not follow from a diagnosis label alone.
Missing or conflicting evidence constrains the accepted research scope.
The mechanism and decision policies remain versioned in `config/track2/`.

## Discovery and public sources

Automated discovery uses admitted public seeds and declared source domains.
Public preparation precedes private execution and receives no controlled case information.
Source records identify query provenance, source versions, and evidence limitations.
Manual literature extraction requires explicit source locators and actual review records.
Software support for manual evidence does not establish that a participant performed such review.

| Source | Declared version or role |
| --- | --- |
| Open Targets Platform | 26.06 target and drug contracts. |
| LINCS | GSE92742, pinned Level 5 matrix and metadata. |
| FDA records | 2026-10-16 snapshot declaration. |
| UNII records and names | 2026-08-04 source files. |
| Public genotype evidence | Optional authenticated export under its declared producer contract. |
| Functional and literature evidence | Versioned functional-source contracts and explicit admitted source references. |

`config/track2/discovery-sources.json` records exact native discovery identities.
The [dependency declaration](../config/installation/dependencies.json) records acquisition and delivery limits.
This table describes supported software inputs.
Submission-specific source use and proprietary sources remain unavailable in this source release.

## Tissue evidence, scoring, and calibration

The exploratory tissue policy requires an accepted loss-of-function mechanism with a declared increase-function objective.
It applies genetic and drug quality rules before comparisons.
The policy selects measured landmark features and preserves cell and tissue provenance.
Paired methods compare compatible overexpression and knockdown evidence.
Single-method evidence retains its separate interpretation.
A failed paired calculation does not automatically select a different method.

Declared strength calculations use retained correlations.
Paired strength is `(r_oe-r_kd)/2`.
Overexpression-only strength is `r_oe`.
Knockdown-only strength is `-r_kd`.
Aggregation combines declared dose, time, cell, and tissue levels under fixed rules.
Contradictory or undefined conditions retain their explicit status.

The derived-scoring policy uses 10000 deterministic background draws with seed 0.
Calibration preserves paired values and declared gene-label randomization.
Percentile is `100*(2*L+E)/(2*B)` under the frozen-background contract.
These values describe relative evidence under that contract.
They do not establish clinical significance or treatment benefit.

Report ranking uses exact percentile fractions and unrounded strength.
Display ties retain stable identities.
Uncalibrated evidence remains in an unranked supplement.
Presentation limits select main tables and figures without deleting complete supplementary evidence.
FDA ingredient eligibility does not establish patient suitability or an accepted clinical use.

## Reports, completion, and recovery

Fresh reporting covers evidence admission, calibration, composition, figures, and finalization.
Complete report artifacts can accompany partial science and retain exit code `2`.
Operational failure cannot establish finalized success.
Actual input, mechanism, and report reviews retain separate identities and scopes.

`track2 report` creates a separate report revision from supported retained run evidence.
Compatible regeneration needs no source VCF, original matrix access, or repeated discovery.
Public scientific resume retains its frozen scientific software and inputs.
Applying later scientific software requires a new public output directory.
A supported separate reporter does not replace that scientific package.

Supported historical transitions require exact producer admission and matching scientific dependencies.
`--recover-unfinished` additionally requires retained production proof.
`--refresh-reporter` applies only before completed source admission or report stage work.
It preserves the previous reporter and failed attempts.
Unknown producers and changed scientific dependencies stop before destructive cleanup.
The [command inventory](portable-command-interface.md) describes these entry points.

## Official workbook coverage

The [dated source record](../release/challenge-sources.json) identifies the official workbook.
Most prompts are voluntary.
Applicable commercial generative-AI disclosure remains required.
These README methods are a project requirement and complement the separate challenge report.

| Workbook topic | Release coverage |
| --- | --- |
| Team | GenoSapiens. |
| Variant-to-medication approach | Mechanism, discovery, tissue evidence, scoring, and report sections above. |
| Automated identification | Declared seed-driven queries and retained perturbation comparisons. |
| Manual review | Source-located manual evidence and actual review records. Submission-specific effort remains unavailable. |
| Public sources | Source table and exact discovery contracts. |
| Proprietary sources | Submission-specific use remains unavailable. |
| Mechanism | Molecular effect, phase, pathways, downstream consequence, and uncertainty above. |
| Time and effort | Shared methods distinguish measured evidence from unavailable runtime, cost, and participant effort. |
| AI disclosure | Exact participant-reported statement in shared methods. |
| Abstract | Abstract above, with strengths and limitations. |

## Glossary

- LINCS: Library of Integrated Network-Based Cellular Signatures.
- FDA: United States Food and Drug Administration.
- UNII: Unique Ingredient Identifier.
- Mechanism: Declared molecular effect and its biological interpretation under retained evidence limits.
- Calibration: Comparison with a declared background distribution.
- Percentile: Relative position under the frozen-background contract.
- `r_oe`: Retained overexpression correlation.
- `r_kd`: Retained knockdown correlation.
- `L`: Number of background strengths below the observed strength.
- `E`: Number of background strengths equal to the observed strength.
- `B`: Number of background draws.
- Producer: Software and contracts that created retained evidence.
- Reporter: Separate software that creates report artifacts from retained evidence.

Paths, identifiers, source names, formulas, and command options retain their exact spelling for traceability.
