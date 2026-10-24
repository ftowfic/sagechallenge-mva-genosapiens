# VCF evidence contracts

D42 introduces a current evidence policy for source parsing and quality decisions.
Shared sample, record, header, and summary evidence binds its SHA-256 digest.
Strict readers reject missing, unknown, or stale policy fields.
Historical reports retain their original policies and producer identities.

The policy requires exact record `FILTER=PASS` and existing sample and variant-class checks for Track 1 quality acceptance.
It separates allele depths from sample depth.
The allele-balance denominator remains FORMAT/DP.
Existing numeric thresholds remain unchanged.
An unexplained AD total above DP requires unresolved quality evidence.

Positive PID grouping requires a documented producer contract bound to the source digest and current phase policy.
The contract defines its version, namespace, block scope, and relationship to GT.
Header declarations alone cannot supply this contract.
Group scope includes input, sample, and contig.
GT/PS and PGT/PID use separate namespaces.

The default header limit is 8 MiB.
The default record and reference-request limits are 4 MiB.
Each source header records its effective read limits.
Inspection limits do not replace these production limits.

Use fresh output for affected construction after a policy change.
Preserve earlier checkpoints and reports.
Reuse unchanged evidence only when authenticated dependencies establish compatibility.

## Terms

- VCF means Variant Call Format.
- QC means quality control.
- GT records genotype alleles.
- AD records allele depths.
- FORMAT/DP records sample depth.
- PS and PID identify phase groups within their respective contracts.
- PGT records a producer-defined phased genotype.
- A namespace separates identifiers with different meanings.
- MiB means 1,048,576 bytes.
- SHA-256 identifies bytes through a cryptographic digest.

Exact field names, policy identifiers, and numeric tokens retain their source forms for traceability.
Structural checks and manual meaning review do not establish full ASD-STE100 dictionary compliance.

## First-read bounds

Header intake and concordance use the same bounded line reader as QC and resumed source parsing.
The reader checks cumulative header bytes before decoding.
It checks each record before splitting fields.
Concordance checks each REF request before reference allocation.
The operator passes its configured limits to concordance.
Concordance reports and source headers retain those limits.
A capacity failure leaves incomplete evidence and produces no successful report.
A failed scan does not read the remaining file merely to complete a digest.

## Sample measurements

The parser maps each sample through its record's ordered FORMAT keys.
It requires GT first when GT exists.
Source evidence distinguishes absent keys, omitted trailing values, explicit dots, and empty tokens.
Malformed consumed values remain errors.
Unused fields retain their exact tokens without typed validation.
This includes PL when no consumer uses its likelihoods.

AD and FORMAT/DP retain separate values and source states.
An AD total above DP does not make VCF syntax invalid.
The default policy marks that relationship unresolved and retains a technical obligation.
It preserves ALT-depth/FORMAT-DP ratios, including values above one.
It does not substitute INFO/DP or summed AD.

Consumed Number=A and Number=R arrays follow ALT cardinality.
Number=G uses combinations from allele count and known ploidy.
Unknown ploidy leaves Number=G cardinality unresolved.
Ploidy means the number of genotype slots.
Cardinality means the expected number of values.
PL records genotype likelihoods in the source format.

## Phase states

Sample evidence distinguishes absent, incomplete, conflicting, unresolved, and usable phase evidence.
A PGT/GT copy conflict remains conflicting even without a producer contract.
A missing PGT or PID partner leaves incomplete evidence.
A complete pair without an applicable producer contract remains unresolved.
Retention does not grant quality acceptance.

Standard phased GT retains explicit PS or its implicit same-contig phase set.
The existing same-record diploid `1/2` rule remains available.
PID group comparisons also require matching contract digests, which include producer versions.
Different source inputs, samples, contigs, or namespaces cannot share a phase group.
Opposing GT/PS and applicable PGT/PID relations remain conflicting.

## Early filtering

The early filter retains phase uncertainty for each affected sample and ALT occurrence.
Protection names the actual PGT/PID or GT evidence and its unresolved reason.
Header declarations alone do not protect a record.
Unresolved AD/DP evidence also retains its technical obligation before common-allele exclusion.
The source ledger retains complete occurrence accounting and original ALT indexes.

## Exact PASS acceptance

Every supplied call requires exact record `FILTER=PASS` for quality acceptance.
A dot retains the distinct unfiltered state.
Named filters retain the failed-filter state.
Sample FT and header declarations cannot supply record PASS.
Existing sample and class-specific checks still apply.
PASS cannot resolve missing measurements or failed metrics.

Admission records each occurrence's quality state, source identifier, sample, and ALT index.
Repeated observations retain separate dispositions.
A PASS occurrence cannot change another occurrence's FILTER value.
Non-PASS calls remain available for source accounting and technical review.
A later read assessment cannot retroactively grant PASS to a supplied occurrence.

## Parsing caches

Each read owner prepares immutable header maps once.
The incremental reducer and source workers also retain their own immutable maps.
Each FORMAT cache retains at most 128 layouts and 256 KiB of charged storage.
The cache uses the complete ordered key tuple.
It checks duplicate keys, GT order, declared keys, and sample arity before interpretation.
Eviction rebuilds a position map without changing evidence.
Unused fields retain raw tokens without typed parsing.

Independent source owners authenticate source bytes again.
A saved stat tuple or supplied digest cannot replace that authentication.
The existing concordance owner closes before QC starts, so it cannot authorize live-owner digest reuse.
Source stability checks remain active throughout each read.
These local changes do not establish a whole-run speed estimate.

## Numeric and reference evidence

Each record retains its original QUAL token, consumed numeric value, raw line, and source digest.
Typed QUAL uses binary64 conversion.
Rendering comparisons identify the adapter version and inspected record count.
A different rendered decimal value does not replace source evidence or establish a genotype or ranking change.
The retained Q25 inspection covers a bounded prefix and two BCFtools versions.
It does not establish whole-file equivalence or behavior for other versions.
Binary64 means a floating-point representation with 64 bits.

Concordance retains every reference declaration in source order.
A declaration without an established operational mapping remains unresolved metadata.
An unresolved path does not prove compatibility or contradiction.
Operational identity still requires the authenticated reference, explicit contig mappings, matching lengths, and complete REF concordance.
The parser does not guess contig aliases or repeat a scan solely to explain a path.
A build label alone cannot authorize a reference join.

## Migration

Current source payloads use version 2.
Current policy bundles use version 3 and include the VCF evidence policy digest.
Source and stage bindings also require that digest.
Admission rows require the same policy field.
Current diagnostic reports use version 5. Current audit manifests use version 6.

Strict readers reject affected predecessor products with missing or stale evidence.
Their files retain their original identities and remain historical evidence.
An old package or checkpoint does not become current through a version-label change.
Unchanged reference and resource artifacts retain their identities when their authenticated dependencies remain valid.

Create fresh affected source, filter, admission, and report outputs under the current policy.
Preserve earlier packages, checkpoints, and reviews.
Do not resume affected scientific construction from a predecessor policy.
Package qualification remains separate from these synthetic implementation checks.

## Engineering package

The unactivated engineering package contains a snapshot of the current shared source tree and its test dependencies.
Its manifest binds each file, the executable inventory, the runtime image, and the VCF evidence policy.
The snapshot includes existing shared-code changes outside D42.
The package does not change the installation selector or establish public qualification.

Producer and consumer checks use tiny independent synthetic fixtures inside the pinned image.
They cover source accounting, filter protection, admission evidence, report publication, and exact cached/uncached parity.
Non-PASS selected support remains unmet even after a synthetic successful read assessment.
The checks also reject stale source, admission, and diagnostic products.
Each invocation has a 900-second limit.

Preserve the package manifest and its test records.
Use the qualified operator workflow before activating a future release.
Do not use this engineering package to resume predecessor scientific outputs.

An executable inventory lists the code and configuration files that define one implementation.
A producer writes evidence. A consumer reads that evidence under its declared contract.
