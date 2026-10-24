# Supported host

The supported host baseline has these components:

- Linux on AMD64.
- Python 3.12 with its standard library.
- Bash and GNU core utilities, including `realpath`, `stat`, and `id`.
- Docker with cgroups v2 and the declared runtime images.

Host launchers need no Python virtual environment.
Image dependencies belong to their image construction contracts.
Other architectures and Python versions have no compatibility claim.

Run commands as an ordinary user with a nonzero user and group identifier.
That user must have permission to run Docker.
The checkout must permit package and installation metadata writes.
The reference root must permit writes to each declared preparation destination.
Runtime source directories must remain readable by that user.

Private files require mode `0400` or `0600`.
Private directories require mode `0700`.
The launcher refuses unsafe symbolic links.
Use paths without line breaks.
Quote paths that contain spaces.

Track 1 uses 12 CPUs and 24 GiB of container memory.
Its memory-plus-swap limit is also 24 GiB.
It requires 26 GiB of available host memory.
It limits execution to 512 process identifiers and a 1 GiB temporary filesystem.

Track 2 uses two CPUs per worker and four GiB of memory by default.
Use `--memory-gib N` for an explicit memory limit.
Its memory-plus-swap limit equals the memory limit.
It limits execution to 256 process identifiers and a 256 MiB temporary filesystem.
Workers execute sequentially under host control.

Both tracks use offline scientific containers with read-only root filesystems.
They remove Linux capabilities and prohibit privilege escalation.
They receive no Docker socket.
Public acquisition remains a separate preparation operation.

The software sets no minimum free-disk reserve.
Preparation needs space for source files, derived resources, and incomplete outputs.
Private execution needs space for checkpoints and report artifacts.
A storage failure preserves incomplete status.
It does not authorize a smaller scientific scope.

## Reference roots and resume

A fresh run defaults to `/mnt/reference`.
Use `--reference-root PATH` to select another installed root.
New runs save the selected root.
Resume restores that root when the option is absent.
An explicit conflicting root fails before preparation or output changes.

Track 1 resume requires preparation contract `track1-early-scope/v1`.
A saved reference root does not make an older preparation contract compatible.
An explicit original root also cannot bypass that rejection.
The launcher rejects the legacy contract before preparation, installation, locks, scratch creation, or run changes.

1. Use a fresh output directory for an incompatible legacy run.
2. Preserve the original installed package and run metadata.
3. Preserve checkpoints and failure evidence for separate recovery.

Compatible preparation contracts use ordinary direct dependency checks.
Preparation, annotation, admission, and science have separate reuse decisions.
Changed scientific dependencies require their actual recomputation scope.
An isolated repair package's science-only restart evidence cannot establish full-release compatibility.
This guide does not establish safe resume of `profile_clinical_query` with `incomplete_required_evidence`.
Admission checkpoint migration through `track1-migrate` remains unsupported.

Older Track 2 runs can recover their root from an authenticated local package.
The saved run must bind that package digest.
Its package manifest must bind `resources.json`.
If that metadata is absent, the command reports the required recovery.

Restore the original package metadata under `.track2-packages` before resuming that older run.
Do not edit saved digests or paths to force acceptance.
A private checkpoint copy on another host has no migration support through these commands.

## Rules reuse

First registration requires an explicit official rules file:

```bash
./track2 setup --rules /private/official-rules.txt \
  --reference-root '/srv/public references'
```

Setup checks missing rules before public preparation.
An existing valid registration retains its authenticated rules copy.
It does not need the original downloaded file.
A conflicting explicit file fails before preparation or registration changes.
Expired registration and stopped incident environments retain their existing recovery requirements.
Setup does not renew consent or reset incident state.

## Glossary

- AMD64: The 64-bit x86 processor architecture.
- Checkout: Directory containing the installed source files.
- Reference root: Host directory containing public resources and their preparation records.
- Digest: Cryptographic value that identifies exact file bytes.
- Manifest: Structured record that binds resource or package files.
- Checkpoint: Completed stage output with authenticated dependencies.
- GiB: 1,073,741,824 bytes.
- MiB: 1,048,576 bytes.
- cgroups v2: Linux controls that limit container resource use.

Paths, identifiers, and option names retain their exact spelling for traceability.
