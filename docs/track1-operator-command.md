# Track 1 operator command

Track 1 prepares its package automatically and runs offline with installed public reference resources.
Ordinary execution and local result publication need no engineering activation, qualification record, or operator override.
External submission remains a separate operator action.

## Installation and private inputs

Use Linux AMD64, Bash, Python 3.12, and Docker with cgroups v2.
Install the public resources before case execution.
The default reference root is `/mnt/reference`.
Use `--reference-root PATH` to select another host directory.
The launcher prepares its package automatically inside Docker.
It uses host Python 3.12 with the standard library only.
You need no host virtual environment, version JSON, package hash, or manual installation pointer.
See [automatic preparation](track1-installation.md) for details.
See the [supported host and recovery guide](portable-host.md) for storage and reference-root requirements.

The first run prepares the public selection index when no compatible completed index exists.
This public-only container receives no case inputs or output directory.
Later compatible runs reuse the completed index and package.
Case processing starts only after preparation succeeds.

Run as the non-root owner of the input files.
Keep input files private with mode `0400` or `0600`.
Use an output directory that you own with mode `0700`.
The launcher creates a missing output directory with that mode.

Keep case inputs, sample names, progress, and results on the local machine.
Do not send these data to an LLM or an external service.
The container has no network access and uses a read-only root filesystem.
It enforces 12 CPUs, 24 GiB memory, no added swap, and 512 process identifiers.
The existing host check requires 26 GiB of available memory.
Use `--allow-low-host-memory` to accept memory exhaustion risk below that threshold for one invocation.
The command records this choice and the observed available memory.
Container limits and scientific dependency checks still apply.
The option does not carry forward to a later resume.
Track 1 has no free-disk reserve requirement.
Actual resource failures stop dependent work and preserve completed checkpoints.

## Start a run

```bash
./track1 run --vcf /private/case.vcf.gz --hpo /private/phenotypes.txt \
  --output /private/results
```

An alternate reference root can contain spaces:

```bash
./track1 run --vcf /private/case.vcf.gz --hpo /private/phenotypes.txt \
  --output /private/results --reference-root '/srv/public references'
```

Preparation and execution use the same selected host directory.
Container resource paths remain unchanged.
See the [command interface inventory](portable-command-interface.md) for supported options.

The command selects the only sample automatically.
For multiple samples, it lists sample names on the private terminal and stops before body processing.
Specify the required sample with `--sample NAME`.
The first sample column has no special status.

Supply the case phenotype terms in the HPO file.
The public HPO ontology supplies term definitions, not patient findings.
Use one identifier per line in UTF-8 text.
Each identifier starts with `HP:` and contains seven digits.

These identifiers show syntax only.

```text
HP:0001250
HP:0001263
```

The parser accepts LF or CRLF line endings and an optional final newline.
It ignores blank lines with a warning and preserves the original bytes and line numbers.
It rejects a file without HPO identifiers.
It also rejects headings, comments, labels, JSON, commas, and free text.
Unsupplied observations remain unknown.

The program checks every source REF against the installed reference before preparation.
It preserves the supplied VCF bytes.
The existing rules for reference declarations, PASS quality, AD/DP evidence, and phase remain active.
See [VCF evidence contracts](track1-vcf-evidence.md) for those rules.
An explicit `--concordance PATH` must identify a complete accepted report for the same input and reference.

Add `--max-hours HOURS` to set a positive finite time limit for the invocation.
The time limit includes automatic public preparation.
Without this option, the invocation has no overall time limit.

## Resume

```bash
./track1 run --output /private/results --resume
```

This package requires the preparation contract `track1-early-scope/v1`.
A run without this contract stops before installation, input processing, locks, or checkpoint cleanup.

1. Use a fresh output directory for a legacy run.
2. Preserve the original package and run directory.
3. Keep the original checkpoints and failure evidence for separate recovery.

This transition does not recover the older `profile_clinical_query` failure with `incomplete_required_evidence`.
Admission checkpoint migration remains unsupported.

For a current-contract run, resume loads omitted inputs from `.operator/run.json`.
An explicit conflicting input fails before output changes.
A nonempty output directory without `--resume` also fails before processing.

The current manifest binds inputs, sample, reference versions, source versions, stage dependencies, and output locations.
Each completed stage records its direct dependencies and output digest.
The program states the restart stage and reason before processing.

A selection change repeats selection and dependent work while it reuses compatible annotation.
A report change reuses compatible scientific outputs.
A different package digest alone does not invalidate every stage.
Missing or invalid stages restart from compatible earlier outputs.
Corrupted completed outputs remain explicit errors.

An output lock prevents concurrent writers.
The worker keeps that lock if its supervisor stops.
The program can replace disposable partial work after it acquires the lock.
Published reports and review records remain intact.

A compatible current-contract resume uses its existing output directory.
It requires no attempt directory, activation receipt, or previous campaign.
`.operator/current` contains the current invocation files.
`.operator/last-error.json` retains the most recent stage error when execution continues.
Old invocation files do not authorize checkpoint reuse.

## Admission diagnostics

Admission checks and records inputs for selection and analysis.
The private terminal shows the failed operation, exception category, rejection reason, and recovery action when current detail exists.
These labels use predefined values.
They exclude exception messages, tracebacks, case identifiers, biological payloads, and paths.

The local record is `analysis/analysis-work/state/admission-diagnostic.json` under the output directory.
It contains private invocation and admission bindings, monotonic times, window counters, and at most four secondary failures.
The maximum record size is 16 KiB.
A terminal diagnostic binds to the current receipt and primary failure.
A missing, stale, or invalid record produces an unavailable-detail notice.
The unsuccessful status remains unchanged.

Keep the diagnostic on the private machine.
Do not share the record, binding digests, or row totals.
A summary for sharing must contain only the predefined failure labels.

| Operation group | Meaning |
| --- | --- |
| `admission_result_collect` | The coordinator collects a worker result. |
| `admission_worker_execute` | A worker checks and encodes admission inputs. |
| `admission_screening_receive` | The coordinator receives a screened window. |
| `profile_window_validate`, `profile_window_decode` | The profile owner checks or decodes a window. |
| `profile_clinical_query`, `profile_clinical_bind`, `profile_clinical_adapt` | The clinical owner queries, checks bindings, or converts evidence. |
| `profile_row_convert`, `profile_row_persist` | The profile owner constructs or saves facts. |
| `admission_window_persist`, `admission_window_retire` | The coordinator saves a window or records its committed state. |
| `admission_complete` | The coordinator checks closure and publishes completion. |

The monitor prints an active-operation notice at most once every 30 seconds.
The notice identifies the last observation and elapsed operation time.
It does not establish forward progress, a failure cause, or a completion estimate.
Elapsed runtime alone does not establish a timeout without an overall deadline.
Existing operation deadlines and resource limits still apply.

Screened, pending, committed, cancelled, and sealed windows describe different states.
For example, 75 screened windows and 74 committed windows leave one pending window.
Cleanup can cancel that window without committing it.
A stopped thread or an empty queue does not establish scientific completion.

| Recovery label | Local action |
| --- | --- |
| `inspect_local_diagnostic_then_resume` | Review the static failure labels before ordinary resume. |
| `check_required_evidence_then_resume` | Check required evidence and source bindings before ordinary resume. |
| `check_storage_then_resume` | Check local storage and private file permissions before ordinary resume. |
| `check_resources_then_resume` | Check runtime resource limits before ordinary resume. |

An incomplete required lookup remains fatal and does not become an empty successful result.
The program does not discard the affected row or retry scientific work automatically.
The original reported failure cause remains unknown.
Counters alone cannot identify that cause.

Use ordinary resume after resolving the local condition in a compatible current-contract run.
Compatible completed preparation and annotation remain reusable after this diagnostic repair.
Incomplete selection and dependent stages restart under their current contracts.
Resume does not promise reuse of committed windows from an unsealed admission product.
Scientific dependency changes still invalidate affected stages.
Diagnostics never authorize completion or treatment decisions.

## Selection and results

Initial quality checks, population decisions, and required annotation still apply to the original input universe.
The shared selector then retains an allele through any applicable route:

- An exact pathogenic, likely pathogenic, or curated damaging functional allele assertion.
- A qualifying existing SIFT or PolyPhen prediction in an established disease gene.
- A supplied structural call within the declared native scope.
- A declared coding or splice consequence on a relevant transcript in an established disease gene.

Track 1 applies this policy across the established public disease-gene catalog without diagnosis input.
Generic regulatory overlap, an unknown effect, or an ordinary synonymous consequence does not independently retain an allele.
Repeat, mitochondrial, and mosaic calls enter the deferred inventory when no route applies.

The compact deferred inventory preserves identities, samples, source occurrences, annotation summaries, reasons, and source versions.
It remains in the ordered profile decision shards.
Deferred alleles do not enter contexts, pairs, scores, mechanisms, or drug inputs.
Retained alleles keep complete relevant pairs, phase uncertainty, partial-member penalties, and selected-member technical checks.

Selection can omit novel noncoding causes and unknown effects without qualifying evidence.
A partner outside the retained universe can still exist in the source genome.
Deferral does not establish benignity or treatment suitability.

The command prints the final report paths.
A local submission CSV appears only when the applicable scientific, reference, technical, and format checks pass.
The program does not submit that file externally.

## Glossary

- Allele: A specific sequence change or a supported native variant identity.
- Source occurrence: One original observation of an allele in the supplied input.
- Retained universe: The alleles that qualify for downstream analysis.
- Deferred inventory: The compact record of alleles without a retention route.
- Checkpoint: A completed stage output with authenticated dependencies.
- Admission: The stage that checks and records inputs for selection and analysis.
- Window: A bounded group of inputs processed together.
- Binding: A digest that identifies the exact invocation, contract, receipt, or failure.
- Monotonic time: A local elapsed-time measurement that does not use calendar time.
- Seal: A record that establishes successful completion under the stage contract.
- Primary failure: The first operational failure observed by its owner.
- Secondary failure: A later failure during reporting or cleanup.
- REF: The reference allele recorded in the VCF.
- HPO: Human Phenotype Ontology.
- SHA-256: The digest algorithm used to identify exact bytes.

Paths, option names, identifiers, and exact format tokens remain unchanged for traceability.
