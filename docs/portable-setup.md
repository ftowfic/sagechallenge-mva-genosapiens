# Public resource setup

The [installation declaration](../config/installation/dependencies.json) is the dependency inventory.
Each entry identifies consumers, role, version, integrity bindings, producer inputs, required records, destinations, and availability.
The [source pins](../config/installation/source-pins.json) bind exact FDA and LINCS source bytes.
These records contain no credentials or controlled inputs.

## Availability limits

The current Track 1 runtime registry manifest was available on 2026-10-19.
This check downloaded no image layers.
The declared Track 2 LINCS, signaturesearch, and discovery-builder image manifests were unavailable at their registry paths.
The signaturesearch base also lacks a complete pinned reconstruction recipe.
Its recorded dependency inventory does not establish reproducible construction.
Fresh Track 2 execution needs an authenticated copy of the exact image or a separately checked replacement binding.

Build recipes and dependency locks appear in each image entry.
A different build digest cannot replace a runtime pin silently.
Follow [the existing image policy](../requirements/supply-chain-policy.json) before binding another digest.
Do not use the broad image publication script for installation.

Several consumer contracts require exact generated manifests from frozen releases.
Their historical builders can produce different timestamps and manifest bytes.
Those dependencies need authenticated frozen artifacts that match the declared digests.
The declaration records unavailable or unverified archives explicitly.
No command invents historical receipts or replacement digests.

The FDA download URL can serve a later snapshot.
Preparation rejects bytes that differ from the frozen October snapshot.
LINCS and Open Targets acquisition code has synthetic fixture coverage.
That coverage does not prove current remote availability.

## Layout and commands

Create a reference root owned by the operator:

```bash
mkdir -p '/srv/public references'
```

The resource declaration gives exact paths beneath that root.
Its main directory roles are:

| Directory | Contents |
| --- | --- |
| `releases` | Frozen reference, annotation, clinical, and ontology sources |
| `manifests` | Authenticated source manifests and local preparation records |
| `derived` | Public indexes, discovery catalogs, and functional resources |
| `downloads` | Public acquisition inputs required by existing builders |
| `lincs/GSE92742` | Frozen metadata, matrix, integrity receipt, and transfer receipt |
| `fda/2026-10-16` | Frozen FDA archives and truthful download metadata |

Inspect a dependency without changing the reference root:

```bash
python3 scripts/prepare_public_resources.py grch38 --describe
```

Prepare metadata without inherited receipt files:

```bash
python3 scripts/prepare_public_resources.py lincs-metadata --reference-root '/srv/public references'
python3 scripts/prepare_public_resources.py fda --reference-root '/srv/public references'
python3 scripts/prepare_public_resources.py open-targets-drugs --reference-root '/srv/public references'
```

Obtain the declared compressed LINCS archive from its source URL under the source terms.
Prepare its matrix and receipts:

```bash
python3 scripts/prepare_public_resources.py lincs-matrix \
  --source /public/downloads/GSE92742-level5.gctx.gz --reference-root '/srv/public references'
```

This explicit operation authenticates and decompresses the archive in one bounded-memory pass.
The filename above illustrates a local archive path.
Its bytes must match the exact compressed and expanded integrity bindings.
Incomplete files remain incomplete after interruption.
Repeat the same command to retry incomplete acquisition.
Compatible completed preparation reuses local completion without another full source hash.

After all required frozen inputs exist, prepare the derived resources:

```bash
python3 scripts/prepare_public_resources.py public-index --reference-root '/srv/public references'
python3 scripts/prepare_public_resources.py discovery --reference-root '/srv/public references'
```

Normal launchers also prepare their current packages automatically.
The explicit commands require no case input, installation pointer, or package cache.
Unavailable frozen roles fail with an installation limit.
The optional public-study export can remain absent during native discovery preparation.

The explicit discovery command authenticates matrix and target inputs before catalog construction.
It reuses valid local matrix preparation without another matrix read.
Normal launchers require local admission or compatible completed catalog evidence.
They report the explicit recovery command when neither exists.

## Relocation and recovery

A relocated checkout derives its paths from the installed launcher.
Its first preparation regenerates local pointers and any affected package metadata.
An alternate reference root supplies all host mount sources consistently.
Existing runs remain bound to their saved reference root.

Copied resources can have different device, inode, and timestamp values.
Copied completion records alone do not authenticate their new local files.
Keep original source manifests and receipts with the copy.
Authenticate copied discovery sources and rebuild affected catalog metadata explicitly:

```bash
python3 scripts/prepare_public_resources.py discovery --admit-copied \
  --reference-root '/srv/copied public references'
```

This operation can read the complete copied matrix once to authenticate its declared digest.
It records new local production separately and preserves original producer evidence.
Later compatible startup reuses the resulting completion.
Never edit stored file identities to authorize copied bytes.
Missing authentication requires restoration of the exact source and its original evidence.

Each new local completion has a separate producer record that binds its bytes and filesystem identity.
Copying or editing the completion invalidates that local binding.
Explicit preparation preserves previous records before it records new local completion.
Existing compatible catalogs retain their original reuse contract.

Other source owners retain their manifest and index checks.
They rebuild affected derived indexes when their existing contracts require it.
These commands do not migrate private checkpoints or change access permissions.

## Glossary

- Digest: Cryptographic value that identifies exact file bytes.
- Receipt: Record of actual acquisition or production and its completion state.
- Manifest: Structured record that identifies resource files and integrity bindings.
- Reference root: Host directory containing public sources and preparation records.
- Local file identity: Device, inode, size, and timestamps recorded on one filesystem.
- Completion: Evidence that the declared production step finished under its bound inputs.

Paths, identifiers, source versions, and structured commands retain their exact spelling for traceability.
