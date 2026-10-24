# sagechallenge-mva-genosapiens

This repository provides two rare-disease research workflows.
Track 1 ranks variants from a VCF and phenotype terms.
Track 2 evaluates genotype mechanisms and medication hypotheses.
These outputs do not establish a diagnosis, treatment, or submission readiness.

Contact: Fadi Towfic ([contact@genosapiens.com](mailto:contact@genosapiens.com)).

Fresh production installation remains limited by unavailable frozen resources and runtime images.
The [setup guide](docs/portable-setup.md) identifies each limit and its recovery requirements.
Synthetic checks establish software behavior only.
They do not establish successful production installation on another host.

## Track 1 methods

Track 1 checks reference compatibility, source accounting, quality, population evidence, and annotation before retained genotype evaluation.
Diagnosis-blind ranking preserves exact variant pairs, incomplete hypotheses, and scientific eligibility limits.
Explicit rapid modes defer contexts without declaring negative findings.
Persistent workers and ordered publication retain dependency-bound checkpoints and shared evidence.
Read the [detailed Track 1 methods](docs/track1-methods.md) for the abstract, sources, workbook coverage, and limits.

## Track 2 methods

Track 2 connects retained genotype mechanisms to medication research hypotheses through declared target, drug, and perturbation evidence.
Tissue and cell comparisons preserve method differences, conflicting evidence, calibration scope, and incomplete outcomes.
Separate reporting software can revise supported retained results while preserving original scientific identities and review records.
Research nominations do not establish treatment benefit or suitability.
Read the [detailed Track 2 methods](docs/track2-methods.md) for the abstract, sources, workbook coverage, and limits.

These README methods complement the separately submitted challenge report.
Comprehensive README methods are a project requirement.
The [dated challenge source record](release/challenge-sources.json) distinguishes official requirements, recommendations, and project choices.

## Setup

Use Linux on AMD64, Python 3.12, Bash, and Docker with cgroups v2.
Host launchers use the Python standard library without a virtual environment.
Follow the [host requirements](docs/portable-host.md) for ownership, permissions, storage, and resource limits.

Read the [dependency declaration](config/installation/dependencies.json) before acquiring resources.
It records exact versions, integrity contracts, producers, destinations, source notices, and availability.
Public preparation receives no case files or credentials.

Request controlled inputs through the [gated dataset page](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data).
Each participant must obtain access and accept the current official rules independently.
Download the official rules through the organizer's access process.
The [source notices](docs/portable-notices.md) cite the dataset access page for this source-only release.
The organizer requires a Synapse citation for later publications using the data.
The organizer's [discussion reply](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/13) explains this scope.

Register Track 2 from your terminal:

```bash
./track2 setup --rules /private/official-rules.txt \
  --reference-root '/srv/public references'
```

Setup requires real participant responses.
It does not grant dataset access.
Valid registration reuses its authenticated rules copy on later setup calls.
An explicit conflicting rules file fails before preparation or registration changes.

## Track 1

```bash
./track1 run --vcf /private/case.vcf.gz --hpo /private/phenotypes.txt \
  --output /private/track1-results --reference-root '/srv/public references'
./track1 run --output /private/track1-results --resume
```

See the [Track 1 operator guide](docs/track1-operator-command.md) for input formats and output paths.
The launcher prepares its package automatically from installed public resources.

## Track 2

Prepare native discovery resources after installing their declared sources:

```bash
./track2 prepare-discovery-resources --reference-root '/srv/public references'
```

An optional frozen public-study export has a separate command:

```bash
./track2 install-public-evidence --source /public/frozen-export \
  --reference-root '/srv/public references'
```

Native discovery does not require that optional export.
Start a diagnosis-guided run:

```bash
./track2 run --vcf /private/case.vcf.gz --diagnosis OMIM:000000 \
  --output /private/track2-results --reference-root '/srv/public references'
./track2 run --output /private/track2-results --resume
```

The diagnosis identifier illustrates syntax only.
It does not identify an accepted diagnosis.
Use `--help` for options, including optional HPO, sample selection, and actual review records.
The [command inventory](docs/portable-command-interface.md) lists every supported host option.

Run a reviewed public genotype through native discovery:

```bash
python3 scripts/run_track2_public.py --bundle /public/genotype-bundle \
  --reviewed-plan /public/genotype-bundle/plan.json --output /public/reconstruction \
  --reference-root '/srv/public references' --prepare-only
python3 scripts/run_track2_public.py --output /public/reconstruction --resume
```

The bundle requires authenticated public sources and actual input and mechanism reviews.
This command does not create human review evidence.
The first command prepares the plan.
The second command performs discovery, rendering, synthesis, overview rendering, and finalization.

The Track 2 report is `OUTPUT/synthesis/controlled-work/report/report.md`.
Its directory includes linked figures, source data, complete results, and supplements.
A report file alone does not establish finalized completion.
Check the final command status and completion records.
Exit code `2` can indicate a partial scientific outcome.
Preserve its stated limits.

Regenerate a report from a supported retained run:

```bash
./track2 report --run /private/existing-run --output /private/new-report
./track2 report --run /private/existing-run --output /private/new-report --resume
```

The source run remains read-only.
Compatible report-only work requires no original VCF, matrix access, or repeated discovery.
`--recover-unfinished` requires supported retained production proof.
`--refresh-reporter` requires `--resume` and failure before completed source admission or report stage work.
It preserves the previous reporter and failed attempts.
Unsupported producer transitions stop without replacing original evidence.
See the [command inventory](docs/portable-command-interface.md) for the entry point and its limits.

## Resume and relocation

Track 1 resume requires `track1-early-scope/v1` and checks each stage's direct dependencies.
A legacy preparation contract requires a fresh output directory and preservation of the original run and package.
This release establishes no safe recovery of the older admission failure or private admission checkpoint migration.
An isolated repair package's compatibility evidence does not transfer to this full release.
Track 2 operator scientific resume follows its direct dependency checks.
Public Track 2 scientific resume retains frozen software and inputs.
Use a new public output directory to apply later scientific software.
Supported reporting transitions can select a separate reporter or create a separate report revision.
Resume restores the saved reference root when the option is absent.
Conflicting explicit reference roots or public inputs fail before output changes.

Checkout relocation, alternate reference paths, and copied public resources have separate contracts.
See [resource recovery](docs/portable-setup.md#relocation-and-recovery).
Private checkpoint migration between hosts remains unsupported.

## Methods and notices

Read the [methods and limits](docs/portable-methods.md) and [challenge notices](docs/portable-notices.md).
Retain [LICENSE](LICENSE) and [ATTRIBUTION.md](ATTRIBUTION.md).
Third-party materials retain their source terms.
Delete controlled data and derivatives from all environments within 30 days after challenge close.
The October 16 announcement specifies November 23, 2026. The [October 13 rules update](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/22726e3eeb1b007c754c01a4d8113d32be913a3c/tabs/rules.py) gives the confirmation address.
Historical rules records and registration bindings retain their original bytes.

Paths, command options, identifiers, and official source titles retain their exact spelling for traceability.
