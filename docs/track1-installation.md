# Run Track 1 with Docker

The normal command prepares its local package automatically.
You do not create version JSON, calculate package hashes, or run a package builder.
The launcher needs Linux on AMD64, Bash, Python 3.12, and Docker with cgroups v2.
Public resources default to `/mnt/reference`.
Use `--reference-root PATH` to select another installed root.
Python uses its standard library only.
You need no host virtual environment.

Run this command from the repository:

```bash
./track1 run --vcf /private/case.vcf.gz --hpo /private/phenotypes.txt \
  --output /private/results
```

Replace the input and output paths with your local paths.
See [the operator guide](track1-operator-command.md) for input permissions, sample selection, and phenotype format.

The launcher performs these steps:

1. Check arguments, local Docker, and installed public directories.
2. Prepare the required public selection index when no compatible completed index exists.
3. Prepare the public early scope when no compatible completed scope exists.
4. Prepare the current package from code, configuration, and public metadata.
5. Save the installation settings automatically.
6. Start offline case processing in Docker.

Public preparation receives no case inputs or result directory.
It reads frozen public sources and writes its public index, early scope, and package cache.
The first index build can take time.
Later compatible runs reuse completed preparation without full-resource scans.
Startup reads the published manifest and a small record that identifies the source versions.
Startup does not hash public data files or repeat index reconciliation.
Routine code edits create a new package and reuse the published index.
New public source versions require a new index.
The case container mounts public resources read-only and has no network access.

The early scope uses the complete accepted Track 1 disease-gene catalog.
Its first build reads existing gene coordinates, linked regulatory intervals, and the completed public selection index.
The shared selection rule determines which exact known disease alleles extend the scope.
The producer requires no additional database download or case input.

The producer writes these assets:

- `<reference-root>/derived/track1-early-scope/<request_sha256>/scope.sqlite`
- `<reference-root>/manifests/track1-early-scope/<request_sha256>.json`
- `<package-store>/early-scope.json`
- `<package-root>/track1-early-scope-pin.json`

The request identifies the policies, source versions, reference, and producer code.
Compatible completed assets require only metadata and indexed binding checks.
Reuse does not certify unread database pages again.
Interrupted production preserves incomplete output and publishes no completed scope pin.

Gene scope includes full Open Targets gene spans and no padding.
These spans can omit transcript regions that VEP recognizes.
Missing gene coordinates remain explicit coverage gaps.
Outside-scope calls and calls with unresolved scope retain source evidence without population queries, annotation, or admission.
Source accounting still covers the complete input.

The launcher stores packages under `.track1-packages` and manages `.track1-installation` itself.
A code change creates a new immutable package when needed.
A failed preparation preserves the previous installation and does not start case processing.
The command requires no engineering activation, qualification campaign, predecessor package, or override file.

Resume a compatible run with this command:

```bash
./track1 run --output /private/results --resume
```

Resume loads saved inputs from the current run manifest.
It reuses compatible stages and preserves accepted reports.
External submission remains a separate operator action.

## Glossary

- Package: Immutable code, configuration, and public version metadata for execution.
- Public selection index: Local clinical and curated evidence for exact allele queries.
- Index reconciliation: Comparison of public projection records with their indexed records.
- Run manifest: Saved input identities, stage dependencies, and output locations.
- Early scope: Public regions and exact alleles that can enter expensive processing.
- Scope pin: A small document that identifies the completed scope manifest and its digest.
- Digest: A value that identifies specific file content.
- Admission: Checks that construct analysis inputs from retained calls and their evidence.

Commands, option names, identifiers, and paths retain their exact form for traceability.
