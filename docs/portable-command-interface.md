# Portable command interface

Run each command from the installed checkout.
An absolute launcher path also works from another working directory.
Quote paths that contain spaces.
Use `--help` to read the parser help.

The table lists supported host options.
Worker options belong to the internal container interface.

| Interface | Options |
| --- | --- |
| `track1 run` | `--vcf`, `--hpo`, `--output`, `--sample`, `--concordance`, `--max-hours`, `--reference-root`, `--allow-low-host-memory`, `--resume`, `--science-mode`, `--science-target-contexts`, `--science-hours`, `--science-report-minutes`, `--science-pilot`, `--science-pilot-minutes`, `--science-workers`, `--science-publication-readers`, `--help` |
| `track2` | `--output`, `--vcf`, `--hpo`, `--concordance`, `--review`, `--mechanism-review`, `--source`, `--rules`, `--diagnosis`, `--sample`, `--resume`, `--require-review`, `--max-hours`, `--memory-gib`, `--lincs-workers`, `--figure-workers`, `--reference-root`, `--settings`, `--help` |
| `track2 report` | `--run`, `--output`, `--resume`, `--recover-unfinished`, `--refresh-reporter`, `--settings`, `--require-review`, `--memory-gib`, `--max-hours`, `--help` |
| `scripts/run_track2_public.py` | `--bundle`, `--reviewed-plan`, `--output`, `--reference-root`, `--memory-gib`, `--lincs-workers`, `--figure-workers`, `--resume`, `--prepare-only`, `--help` |
| `scripts/prepare_public_resources.py RESOURCE` | `--reference-root`, `--source`, `--admit-copied`, `--describe`, `--memory-gib`, `--help` |

The Track 2 host parser supports these commands:

- `setup` accepts an official rules file for first registration.
- `run` accepts case inputs and a diagnosis for a new run.
- `prepare-discovery-resources` prepares public resources without case inputs.
- `install-public-evidence` accepts the declared optional export through `--source`.
- `report` creates a separate report revision from supported retained run evidence.

The reference root defaults to `/mnt/reference` for a fresh installation.
It contains the declared public resources and their preparation records.
Changing this path does not authorize private checkpoint migration.
The explicit preparation command lists resource names in its help.
Use `--describe` to inspect a dependency without preparation.
Unavailable frozen dependencies remain installation limits.

## Track 1 science options

An explicit `--science-mode` activates persistent science execution options.
Supported modes are `exhaustive`, `rapid_phenotype`, and `rapid_strict_association`.
Default execution without this option uses a separate path.
Worker and publication-reader options alone do not activate the explicit science path.

`--science-workers` and `--science-publication-readers` each accept 1, 2, 4, or 8.
Each option defaults to one process.
These controls do not change scientific scope or scoring.
The retained four-worker measurement does not change these defaults.
Eight-process options have no completed eight-process benchmark claim.

`--science-target-contexts` defaults to 1500 for rapid scope.
Protected contexts and boundary ties can exceed that target.
Rapid deferral means unevaluated work.
Rapid completion does not establish exhaustive completion.

`--max-hours` sets an invocation deadline that includes public preparation.
Omission of that option sets no overall invocation deadline.
Explicit science options default to a five-hour cumulative budget with a thirty-minute report reserve.
`--science-pilot` selects a twenty-minute pilot.
`--science-pilot-minutes 40` requests the supported explicit extension.
Persistent budget records retain their original time basis across resume.

The internal operator supports `--science-no-budget` and `--science-fresh-trial`.
The standard `track1` host launcher rejects these options.
Elapsed-only timing and fresh-trial timing therefore require their supported internal or development entry points.
They do not become standard host commands through this guide.
Historical validation launchers, specialized repair launchers, and benchmark scripts remain development tools.

## Track 2 report options

`track2 report --run RUN --output REVISION` uses its own parser and dispatch.
It selects reporting software separately from frozen scientific software.
The source run receives a read-only mount.
Compatible reporting requires no source VCF, original matrix access, or repeated discovery.
`--resume` reuses compatible completed report work.

`--recover-unfinished` requires supported historical producer and production evidence.
It preserves missing original completion records and creates separate completion evidence for the revision.
`--refresh-reporter` requires `--resume` and failure before completed source admission or report stage work.
It preserves the previous reporter and failed attempts.
It cannot replace completed report work.
Unsupported producers or changed scientific dependencies stop the transition.

The separate report and its sibling `.operation` directory retain revision state.
Keep both directories for compatible retry.
Original packages, checkpoints, reports, and actual reviews remain preserved.
Scientific incompleteness retains exit code `2` after report completion.
Report review remains separate from original input and mechanism reviews.

## Examples

Register with the downloaded official rules:

```bash
./track2 setup --rules /private/official-rules.txt \
  --reference-root '/srv/public references'
```

Prepare the discovery resources:

```bash
./track2 prepare-discovery-resources --reference-root '/srv/public references'
```

Select explicit rapid science with four workers and four publication readers:

```bash
./track1 run --vcf /private/case.vcf.gz --hpo /private/phenotypes.txt \
  --output /private/rapid-results --reference-root '/srv/public references' \
  --science-mode rapid_phenotype --science-workers 4 --science-publication-readers 4
```

Create and resume a separate report revision:

```bash
./track2 report --run /private/existing-run --output /private/new-report
./track2 report --run /private/existing-run --output /private/new-report --resume
```

Install the optional declared public export:

```bash
./track2 install-public-evidence --source /public/frozen-export \
  --reference-root '/srv/public references'
```

Start a diagnosis-guided run:

```bash
./track2 run --vcf /private/case.vcf.gz --diagnosis OMIM:000000 \
  --output /private/track2-results --reference-root '/srv/public references'
```

The diagnosis identifier illustrates syntax only.
It does not identify an accepted diagnosis.

## Glossary

- Checkout: Directory containing the installed source files.
- Reference root: Host directory containing public resources and their preparation records.
- Checkpoint: Completed stage output with authenticated dependencies.
- Export: Frozen public evidence directory accepted by its declared integrity contract.

Paths, identifiers, and option names retain their exact spelling for traceability.
