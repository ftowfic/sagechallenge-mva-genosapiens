# Methods and limits

Track 1 normalizes and annotates variants under frozen reference and execution contracts.
Its ranking combines declared molecular and phenotype evidence.
Its outputs retain scientific eligibility and coverage limits.
Portability changes do not change the v5 execution limits or scientific policies.

Track 2 accepts diagnosis-guided VCF inputs or reviewed published genotypes.
It retains supported alleles before genotype and mechanism construction.
It separates deferred alleles, uncertain mechanisms, research nominations, and submission readiness.
Optional HPO terms supply context.
They do not replace the allele retention rules.

Native discovery uses frozen public sources and cell-specific LINCS comparisons.
It preserves failed jobs, conflicting conditions, and incomplete evidence.
Report scores describe retained evidence for review.
They do not establish efficacy, safe exposure, or treatment suitability.
Complete exports and supplements preserve evidence beyond the main tables.

Scientific containers have no network access.
Public acquisition precedes case processing and receives no controlled case information.
Source metadata and package fingerprints identify the exact inputs used.
Compatible resume reuses completed work through direct dependency checks.
Finalization checks report artifacts and preserves actual scientific completion status.

Participant-reported AI use:

1. VCF structure assessment: OpenAI API, model identifier Luna/6, commercial terms. Data sharing disabled. No training on inputs or outputs. No feedback provided.

2. Coding: OpenAI Pro 200 plan. The coding assistant had no access to input, output, or controlled data. Data sharing disabled. Improve model unchecked. No feedback provided.

Actual human reviews retain their recorded scope.
This document claims no new scientific review or patient analysis.

## Reproducibility and sources

`RELEASE.json` binds both tracks to one source commit and an exact exported file inventory.
The [Track 1 methods](track1-methods.md) describe preprocessing, annotation, ranking, pairs, and scientific scope.
The [Track 2 methods](track2-methods.md) describe mechanisms, discovery, tissue evidence, scoring, and calibration.
The [command inventory](portable-command-interface.md) identifies supported entry points and option defaults.
The [dependency declaration](../config/installation/dependencies.json) records versions, integrity contracts, producers, and installation limits.

Declared public source support does not establish actual submission use.
Submission-specific proprietary source use, downstream manual curation, and participant effort remain unavailable in this source release.
Actual input, mechanism, and report review records retain their separate scope.
Source publication creates no human scientific acceptance or challenge submission.

## Runtime, cost, and effort evidence

The following retained benchmark used twelve synthetic contexts with production resources.
It used one synthetic TRIP13 variant and a fixed 64-term phenotype profile.
Those terms do not describe one patient.
Every row used package `f2373442f19c61a33103295251942756536b2f8d43467c45d0b69686e24bc458`.
The benchmark operated under twelve-CPU and 24-GiB limits.

| Context workers | Publication readers | Complete seconds | Context seconds | Publication seconds | Sampled peak GiB |
| --- | --- | ---: | ---: | ---: | ---: |
| 1 | 1 | 1448.545 | 269.162 | 332.966 | 16.129 |
| 2 | 2 | 1289.935 | 214.749 | 209.139 | 14.339 |
| 4 | 4 | 1205.765 | 144.044 | 206.070 | 17.924 |

The four-worker configuration completed 84.170 seconds earlier than the two-worker configuration.
Its complete-run speedup over one worker was 1.201×.
This result remains below the 7–10× target.
The benchmark did not control filesystem cache state.
Sampled counters can omit the final sampling interval.
Eight workers and eight readers remain unmeasured under the later user instruction.
Installed process-count options still default to one.

These measurements establish neither whole-genome throughput nor compatibility for a different package.
The retained benchmark report identifies its original package and measurement scope.
Release evidence must separately bind the selected export and tested runtime.
New source publication requires no repetition of these production measurements.

Production scientific comparison and independent report replay remain omitted under the 2026-10-24 authorization.
Omitted checks do not count as passed checks.
The final WGS run did not start under the latest stop instruction.
A future user request is necessary for that continuation.

No matched whole-genome runtime measurement supports this release.
Track 2 production runtime, monetary cost, and participant effort remain unavailable in this source record.
The document supplies no invented estimates for those quantities.

## Glossary

- VCF: Variant Call Format file containing observed sequence variants.
- HPO: Human Phenotype Ontology terms describing observed features.
- LINCS: Library of Integrated Network-Based Cellular Signatures.
- Genotype: Observed allele combination under its measurement and phase limits.
- Fingerprint: Digest of the dependencies that control package or stage reuse.
- Finalization: Checks and publication that establish terminal report status.
- Context: Scientific evaluation group with complete retained candidate membership.
- Publication reader: Process that prepares saved context evidence for ordered report publication.
- GiB: Gibibyte, equal to 1073741824 bytes.

Identifiers and official source titles retain their exact spelling for traceability.
