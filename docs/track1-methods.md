# Track 1 methods

This document describes released software capabilities.
It does not report a patient result or certify an actual submission.
`RELEASE.json` identifies the selected source and exported files.
The [shared methods](portable-methods.md) identify disclosure, measurement, and review limits.

## Abstract

Track 1 ranks supplied sequence variants with declared molecular, population, clinical, and phenotype evidence.
The workflow checks reference compatibility, quality, normalized identity, and source occurrence accounting.
Public preparation precedes private execution.
Early population and scope rules reduce annotation work under explicit protection and uncertainty rules.
VEP annotation and exact public evidence support allele retention.
The workflow constructs genotype hypotheses from complete retained context membership.
It preserves compound heterozygous pairs and identifies incomplete recessive hypotheses separately.
Primary ranking remains diagnosis blind.
EPCR supplies an uncalibrated relative ordering for eligible predictions.

Explicit science modes select exhaustive or rapid context evaluation.
Rapid selection preserves molecular protection, uncertainty, boundary ties, and complete context membership.
Deferred contexts remain unevaluated.
Persistent workers evaluate contexts, while ordered publication combines their retained evidence.
Shared resource bundles support saved checkpoint reading under ordinary dependency checks.
Failed required evidence prevents complete scientific publication.

Source identities, explicit uncertainty, and deterministic evidence support reproducibility.
Incomplete public knowledge, retained scope limits, and missing family information constrain interpretation.
Synthetic execution checks establish software behavior within their declared scope.
They do not establish diagnostic sensitivity, treatment suitability, or whole-genome performance.

## Inputs and preprocessing

Inputs include a VCF, HPO identifiers, a selected sample, and installed public resource contracts.
An optional concordance record must bind the same source and reference.
The workflow preserves supplied VCF bytes and original occurrence mappings.
It checks source REF alleles against the reconstructed GRCh38 reference.
Normalization identifies duplicate alleles without discarding their source occurrences.

Quality, genotype, depth, allelic-depth, phase, and representation rules retain their declared limits.
Missing family, onset, and phase information remains unknown.
The early scope contract binds population, clinical protection, and exception evidence.
Unknown evidence cannot become a complete negative result.
The [VCF evidence guide](track1-vcf-evidence.md) describes these input rules.

## Annotation and evidence

VEP 116 provides transcript consequences and applicable molecular predictions under pinned cache and reference identities.
The focused retention policy combines all applicable routes.
These routes include exact clinical evidence, curated functional evidence, eligible molecular consequences, and supported native structural calls.
An allele can qualify through any applicable route.
Diagnosis labels and phenotype similarity do not select Track 1 alleles.
Generic regulatory overlap and unknown effects do not qualify alone.

Deferred alleles retain reasons, occurrence accounting, and policy identities.
They enter no genotype contexts or scores without qualifying retention evidence.
The focused profile can omit novel noncoding causes and unknown effects.
Required resource failures stop dependent work.
Optional unavailable evidence retains that status.

| Public source | Declared release or contract |
| --- | --- |
| GRCh38 reference | `config/reference-identities/challenge-grch38-masked-v2-hs38d1.json` |
| VEP | 116, `config/public-resources/ensembl-vep-116-grch38.json` |
| gnomAD | 4.1, joint projection and all-sites contracts |
| ClinVar | GRCh38, 2026-08-22 |
| Open Targets Platform | 26.06 |
| Rare-disease evidence | `config/public-resources/rare-disease-evidence-2026-08-30.json` |
| Exomiser resources | 15.1.0, data release 2512 |
| HPO and other ontologies | `config/public-resources/core-ontologies.json` |

The [dependency declaration](../config/installation/dependencies.json) records exact integrity bindings and delivery limits.
Source support does not establish actual use in a submission.
Submission-specific proprietary source use remains unavailable in this source release.

## Contexts, ranking, and pairs

Contexts group a disease, target, requirement, and complete retained candidate membership.
The workflow evaluates genotype hypotheses and records their evidence and eligibility.
EPCR orders eligible predictions without a calibrated causal probability claim.
The competitive representation contains at most ten primary singleton or unordered-pair rows.

Complete compound heterozygous pairs occupy one row when both alleles pass representation and eligibility checks.
An eligible incomplete recessive hypothesis can supply one explicitly labeled partial-member prediction.
Its unresolved partner, incompleteness penalty, and selected-member checks remain explicit.
A partial-member prediction does not establish a complete genotype.

Competitive output remains primary-only.
Other assessed candidates remain in retained evidence and ordering records.
Incidental-finding return and public disclosure require separate resolution.
This software description establishes no family-contact authorization or consent for public case disclosure.

## Science execution and recovery

Default execution and explicit `--science-mode exhaustive` use separate implementation paths.
An explicit science mode activates the persistent science options.
Worker and publication-reader options alone do not activate that path.
Host examples and timer limits appear in the [command inventory](portable-command-interface.md).

`rapid_phenotype` prioritizes reached contexts with phenotype evidence and protected molecular evidence.
`rapid_strict_association` adds its declared direct, strong association requirement.
The policy preserves explicit deferral reasons and uncertainty.
The context target defaults to 1500 and does not split context membership.
Protected contexts and exact boundary ties can exceed that target.
Rapid completion does not establish exhaustive completion.

Workers keep private query sessions and evidence writers.
The coordinator consumes completions in original context order.
Publication readers prepare saved projections for one ordered writer.
Shared checkpoints retain immutable resource references for saved reading after the producing process closes.
Recovery rejects missing, changed, foreign, or incomplete required bundles.
Changed scientific dependencies require their actual recomputation scope.

Preparation, annotation, admission, and science have separate reuse decisions.
An isolated repair package's recovery evidence cannot establish compatibility for the full release.
Preserve original packages, metadata, checkpoints, and failure evidence before any necessary recomputation.
Admission checkpoint migration remains unsupported.
The [host guide](portable-host.md) states the legacy preparation rejection rule.

## Official workbook coverage

The [dated source record](../release/challenge-sources.json) identifies the official workbook.
Most workbook prompts are voluntary.
The commercial generative-AI disclosure is conditional and required when applicable.
Complete README methods are a project requirement.
The README and these methods complement the separately submitted challenge report.

| Workbook topic | Release coverage |
| --- | --- |
| Team and approach | GenoSapiens, Track 1 software capabilities. Submission model numbering remains submission-specific. |
| Detailed approach | Preprocessing, annotation, contexts, ranking, and execution sections above. |
| Automated output | Software generates evidence and prediction artifacts. Submission-specific downstream curation remains unavailable. |
| Manual review | Actual review records retain their scope. This source release establishes no new patient review. |
| Public data | Source table and integrity-bound dependency declaration. |
| Proprietary data | Submission-specific use remains unavailable. |
| Variant pairs | Complete pairs and explicitly incomplete partial-member predictions above. |
| Incidental findings | Primary competitive output and separate disclosure limits above. |
| Runtime and cost | Shared methods separate retained measurements from unavailable production and cost evidence. |
| AI disclosure | Exact participant-reported statement in shared methods. |
| Abstract | Abstract above, with strengths and limitations. |

## Glossary

- VCF: Variant Call Format file containing supplied sequence calls.
- HPO: Human Phenotype Ontology identifiers for observed features.
- VEP: Ensembl Variant Effect Predictor.
- EPCR: Uncalibrated relative value that orders eligible predictions.
- Context: Scientific evaluation group with complete retained member identities.
- Admission: Checks that determine which evidence or candidates can enter a stage.
- Compound heterozygous pair: Two distinct alleles considered together in a recessive genotype hypothesis.
- Checkpoint: Durable completed work with its dependency bindings.

Paths, identifiers, source names, and command options retain their exact spelling for traceability.
