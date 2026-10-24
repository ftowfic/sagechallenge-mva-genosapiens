"""Offline materialization of public model-state validation perturbations.

Selection, calls, phenotype terms and truth stay in deterministic tooling.
Only a static completion status is emitted. Generated calls are synthetic
research probes, not assertions that a public allele causes a disease.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import socket
import subprocess
from types import SimpleNamespace

from genosapiens.grch38_reference import IndexedFasta
from genosapiens.reference_concordance import load_reference_identity
from genosapiens.track1_roster import REQUIRED_MODEL_COVERAGE


REVIEW_STATUSES = (
    "criteria provided, multiple submitters, no conflicts",
    "practice guideline", "reviewed by expert panel",
)
ROSTER_PATH = Path("tests/fixtures/track1/roster/roster.json")
DERIVATION_PATH = Path("config/track1/public-population-development-derivation.json")


def _bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def _write(path, content):
    if path.is_symlink() or any(parent.is_symlink() for parent in path.parents):
        raise ValueError
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)


def _binding(root, path, role=None):
    content = path.read_bytes()
    value = {"path": path.relative_to(root).as_posix(), "size_bytes": len(content), "sha256": _sha(content)}
    if role is not None:
        value["role"] = role
    return value


def _query(connection, sql, parameters):
    return connection.execute(sql, parameters).fetchall()


def _context_choices(connection, references):
    targets = str(references / "releases/open-targets/26.06/output/target/*.parquet")
    eva = str(references / "releases/open-targets/26.06/output/evidence_eva/*.parquet")
    rare = references / "derived/rare-disease-evidence/orphadata-2026-06-23_gencc-2026-08-30_ot-26.06_v6"
    return _query(connection, '''
        SELECT DISTINCT t.id, t.approvedSymbol, t.genomicLocation.chromosome,
               t.genomicLocation.start, t.genomicLocation.end, g.source_disease_id
        FROM read_parquet(?) t JOIN read_parquet(?) g ON g.ensembl_gene_id=t.id
        JOIN read_parquet(?) p ON p.source_disease_id=g.source_disease_id
        WHERE t.genomicLocation.chromosome IN ('1','2','3','4','5','6','7','8','9','10',
             '11','12','13','14','15','16','17','18','19','20','21','22')
          AND t.genomicLocation.end-t.genomicLocation.start BETWEEN 100 AND 500000
          AND p.canonical_hpo_id IS NOT NULL AND coalesce(p.qualifier_not,false)=false
          AND t.id IN (SELECT targetId FROM read_parquet(?)
             WHERE list_contains(allelicRequirements,'Autosomal recessive inheritance')
             AND confidence IN (?,?,?))
        ORDER BY t.id, g.source_disease_id
        ''', [targets, str(rare / "gene_assertion.parquet"), str(rare / "phenotype_assertion.parquet"), eva, *REVIEW_STATUSES])


def _variants(references, chromosome, start, end):
    bcf = references / "derived/gnomad/4.1/joint-frequency/gnomad-joint-v4.1-frequency-projection-v1" / f"gnomad.joint.v4.1.chr{chromosome}.bcf"
    result = subprocess.run(("bcftools", "query", "--regions", f"chr{chromosome}:{start}-{end}",
        "--format", "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/GN_FILTERS_DEFINED\\t%INFO/GN_FILTERS\\t%INFO/GN_AC\\t%INFO/GN_AN\\t%INFO/GN_NHOMALT\\t%INFO/GN_NHEMI\\n", str(bcf)),
        capture_output=True, timeout=300)
    if result.returncode or len(result.stdout) > 64 * 1024**2:
        raise ValueError
    rows = []
    for line in result.stdout.decode("ascii").splitlines():
        chrom, pos, ref, alt, filters_defined, filters, ac, an, hom, hemi = line.split("\t")
        if chrom != "chr" + chromosome:
            raise ValueError
        # BCF transport PASS is synthetic. Only explicit genome source filters
        # and consistent genome counts can qualify a new development probe.
        if (len(ref) == len(alt) == 1 and ref in "ACGT" and alt in "ACGT" and ref != alt
                and filters_defined == "1" and filters in {".", "PASS"}
                and ac.isdecimal() and an.isdecimal() and 0 < int(ac) < int(an)
                and all(count == "." or count.isdecimal() for count in (hom, hemi))
                and 2 * (0 if hom == "." else int(hom)) + (0 if hemi == "." else int(hemi)) <= int(ac)):
            rows.append((chromosome, int(pos), ref, alt, int(ac), int(an)))
    return rows


def _unasserted(connection, references, chromosome, rows):
    # Exclude any allele with an Open Targets EVA record, regardless of strength
    # or condition linkage. No clinical absence is turned into a benign label.
    eva = str(references / "releases/open-targets/26.06/output/evidence_eva/*.parquet")
    asserted = {r[0] for r in _query(connection, 'SELECT DISTINCT variantId FROM read_parquet(?) WHERE starts_with(variantId,?)', [eva, chromosome + "_"])}
    return [row for row in rows if "_".join(map(str,row[:4])) not in asserted]


def _select(connection, references, choices, forbidden_genes, *, dev_chromosome=None):
    for gene_id, symbol, chromosome, start, end, disease in choices:
        if symbol in forbidden_genes or chromosome == dev_chromosome:
            continue
        rows = _unasserted(connection, references, chromosome, _variants(references, chromosome, start, end))
        rows.sort(key=lambda row: (Fraction(row[4], row[5]), row[:4]))
        if len(rows) < 3 or Fraction(rows[1][4], rows[1][5]) >= Fraction(rows[-1][4], rows[-1][5]):
            continue
        rare = references / "derived/rare-disease-evidence/orphadata-2026-06-23_gencc-2026-08-30_ot-26.06_v6"
        hpo = tuple(r[0] for r in _query(connection, '''SELECT DISTINCT canonical_hpo_id FROM read_parquet(?)
            WHERE source_disease_id=? AND canonical_hpo_id IS NOT NULL AND coalesce(qualifier_not,false)=false
            ORDER BY canonical_hpo_id''', [str(rare / "phenotype_assertion.parquet"), disease]))
        if hpo:
            return (gene_id, symbol, chromosome, disease, hpo, (rows[0], rows[1], rows[-1]))
    raise ValueError


def _development_derivation(repository, references, selected):
    low = max(Fraction(row[4],row[5]) for row in selected[-1][:2])
    high = Fraction(selected[-1][2][4],selected[-1][2][5])
    # Exact rational separator in log-frequency space, rounded upward at a
    # declared decimal precision. Neither held-out frequencies nor performance
    # enters this fixed choice. This is stress-test policy selection, not a
    # clinical maximum credible allele frequency calculation.
    from decimal import Decimal, localcontext
    with localcontext() as context:
        context.prec=60
        midpoint=(Decimal(low.numerator)/Decimal(low.denominator)*Decimal(high.numerator)/Decimal(high.denominator)).sqrt()
        scale=10**12
        numerator=int(midpoint*scale)+1
    threshold=Fraction(numerator,scale)
    if not low < threshold < high:
        raise ValueError
    from genosapiens.gnomad_publication import load_gnomad_publication_contract
    from genosapiens.track1_public_validation import (
        GNOMAD_JOINT_MANIFEST_SHA256, OPEN_TARGETS_MANIFEST_SHA256,
        RARE_DISEASE_BUNDLE_SHA256, _bound_bytes,
    )
    publication = load_gnomad_publication_contract(repository / "config/public-resources/gnomad-joint-publication.json")
    manifests = (
        (publication.manifest_path, GNOMAD_JOINT_MANIFEST_SHA256),
        ("manifests/open-targets/open-targets-platform-26.06.0cae4540cfea.json", OPEN_TARGETS_MANIFEST_SHA256),
        ("manifests/rare-disease-evidence/rare-disease-evidence-2026-08-30-v6.cded6e513755.bundle.json", RARE_DISEASE_BUNDLE_SHA256),
    )
    return {
        "schema":"genosapiens.track1-public-population-development-derivation",
        "classification":"public", "validation_only":True, "controlled_use_authorized":False,
        "selection_authority":"development-source-only; no held-out frequencies or evaluator outcomes",
        "algorithm":"geometric midpoint of the two selected development low-frequency upper bound and selected high-frequency AF; upward rational rounding at 12 decimal places",
        "source_manifest_sha256s":sorted(_sha(_bound_bytes(references / path, digest)) for path, digest in manifests),
        "development_selected_low_frequency_upper_bound":str(low), "development_selected_high_frequency":str(high),
        "selected_frequency_threshold":str(threshold),
        "minimum_observed_development_an":min(row[5] for row in selected[-1]),
        "selected_subset":"genomes", "selected_stratum":"adjusted", "selected_ancestry":None,
        "selected_low_identity_sha256s":sorted(_sha(':'.join(map(str,row[:4])).encode('ascii')) for row in selected[-1][:2]),
        "selected_high_identity_sha256":_sha(':'.join(map(str,selected[-1][2][:4])).encode('ascii')),
        "limitations":["Source-selected synthetic stress probes are not independent patients or clinical calibration.", "No curated exception or clinical causal assertion is inferred from selection."],
    }


def materialize_model_state_roster(repository, references):
    if any(name != "lo" for _,name in socket.if_nameindex()):
        raise ValueError
    import duckdb
    connection=duckdb.connect(":memory:")
    connection.execute("SET threads=2")
    connection.execute("SET memory_limit='2GB'")
    roster_path=repository/ROSTER_PATH
    roster=json.loads(roster_path.read_bytes())
    if any(case['generator']['version'] != '1' for case in roster['cases']):
        raise ValueError  # This command never regenerates/relabels its own outputs.
    original_roster_sha=_sha(roster_path.read_bytes())
    choices=_context_choices(connection,references)
    genes={split:{gene for case in roster['cases'] if case['split']==split for entry in case['annotation_expectations'] for gene in entry.get('gene_symbols',())} for split in roster['splits']}
    dev=_select(connection,references,choices,genes['held_out'])
    # Derive before even selecting the fresh held-out source probes.
    derivation=_development_derivation(repository,references,dev)
    derivation_bytes=_bytes(derivation)
    _write(repository/DERIVATION_PATH,derivation_bytes)
    held=_select(connection,references,choices,genes['development']|{dev[1]},dev_chromosome=dev[2])
    contract=load_reference_identity(repository/'config/reference-identities/challenge-grch38-masked-v2-hs38d1.json')
    role=contract.reference_role('public_annotation_benchmark')
    reference=IndexedFasta(references/contract.expected_output(role.fasta_output_role).path,references/contract.expected_output(role.fai_output_role).path)
    from genosapiens.track1_public_validation import (
        _authenticated_public_reference_aliases,
    )
    manifest_path = references / (
        'manifests/grch38/'
        'challenge-grch38-masked-v2-plus-hs38d1.a92c50f61c1b.json'
    )
    source_aliases = _authenticated_public_reference_aliases(
        identity_path=repository/'config/reference-identities/challenge-grch38-masked-v2-hs38d1.json',
        reference_manifest_path=manifest_path,
        reference_root=references,
        reference_contract=contract,
        reference=reference,
    )
    reference_order = {name: index for index, name in enumerate(reference.index)}
    truth_path=repository/roster['causal_truth']['path'];truth=json.loads(truth_path.read_bytes())
    annotations_path=repository/roster['annotation_expectations']['path'];annotations=json.loads(annotations_path.read_bytes())
    roster['required_model_coverage']=list(REQUIRED_MODEL_COVERAGE)
    for case in roster['cases']: case['required_model_coverage']=[]
    coverage=(
        ('unknown-pair', (0,1), False, ['clinical-unasserted-exploratory-pair','completed-population-lookup-unknown-rescue','unknown-count-hypothesis']),
        ('excluded-partner', (0,2), False, ['completed-population-lookup-unknown-rescue','population-excluded-unresolved-partner']),
        ('cis-contradiction', (0,2), True, ['independent-non-population-contradiction']),
    )
    for split,selected in (('development',dev),('held_out',held)):
        gene_id,symbol,chromosome,disease,hpo,rows=selected
        group_id='public-model-gene-group-'+_sha((gene_id+'\n').encode())[:16]
        for case in roster['cases']:
            if case['split']==split and any(symbol in entry.get('gene_symbols',()) for entry in case['annotation_expectations']):case['group_id']=group_id
        baseline=next(case for case in roster['cases'] if case['split']==split and case['scope']=='compound_heterozygous')
        baseline_vcf=repository/next(item['path'] for item in baseline['runner_inputs'] if item['role']=='small_variant_vcf')
        baseline_lines=baseline_vcf.read_text(encoding='ascii').splitlines()
        background_set={'1:602113:T:TGCCCA','1:604358:G:C','1:604360:T:C','1:631859:CG:C','1:1735008:G:GGTTT','1:1735008:G:GTTTT'}
        header=[line for line in baseline_lines if line.startswith('#')]
        # The new calls may use a contig absent from the old scope fixture's
        # valid subset dictionary. Author a new subset dictionary for these new
        # synthetic bytes, preserving existing declarations and public lengths.
        contig_lines = [line for line in header if line.startswith('##contig=')]
        contig_names = {re.search(r'ID=([^,>]+)', line).group(1) for line in contig_lines}
        if chromosome not in contig_names:
            contig_lines.append('##contig=<ID=' + chromosome + ',length=' + str(reference.index[source_aliases[chromosome]].length) + '>')
        contig_lines.sort(key=lambda line: reference_order[source_aliases[re.search(r'ID=([^,>]+)', line).group(1)]])
        header = [line for line in header if not line.startswith('##contig=')]
        insertion = next(index for index, line in enumerate(header) if line.startswith('##reference=')) + 1
        header[insertion:insertion] = contig_lines
        background=[line for line in baseline_lines if not line.startswith('#') and any(':'.join((line.split('\t')[0],line.split('\t')[1],line.split('\t')[3],alt)) in background_set for alt in line.split('\t')[4].split(','))]
        if not any(line.startswith('##FORMAT=<ID=PS,') for line in header):header.insert(-1,'##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Synthetic phase set">')
        for number,(suffix,indices,phased,required) in enumerate(coverage,1):
            case_id=split.replace('_','-')+'-model-state-'+str(number).zfill(2)
            case=deepcopy(baseline);case.update(case_id=case_id,group_id=group_id,truth_id=case_id+'-truth',required_model_coverage=sorted(required),generator={'id':'track1-model-roster-generator','version':'2','seed':101 if split=='development' else 303})
            vcf_path=roster_path.parent/(case_id+'.vcf'); hpo_path=roster_path.parent/(case_id+'.hpo.txt')
            calls=list(background);members=[];expectations=[]
            for index in indices:
                chrom,pos,ref,alt,ac,an=rows[index]
                if reference.fetch('chr'+chrom,pos,len(ref))!=ref:raise ValueError
                member=':'.join(map(str,(chrom,pos,ref,alt)));members.append(member)
                calls.append('\t'.join(map(str,(chrom,pos,'.',ref,alt,'.','PASS','.', 'GT:PS' if phased else 'GT','0|1:1' if phased else '0/1'))))
                expectations.append({'member':member,'gene_symbols':[symbol],'most_severe_consequence':None,'split':split,'assessment_status':'public-target-interval-expectation; actual-VEP-verification-required'})
            calls.sort(key=lambda line:(int(line.split('\t')[0]),int(line.split('\t')[1]),line))
            _write(vcf_path,('\n'.join(header+calls)+'\n').encode('ascii'));_write(hpo_path,('\n'.join(hpo)+'\n').encode('ascii'))
            case['runner_inputs']=[_binding(repository,hpo_path,'hpo_text'),_binding(repository,vcf_path,'small_variant_vcf')]
            format_path=repository/baseline['source_format_bindings'][0];source_format=json.loads(format_path.read_bytes());source_format['input_sha256']=_sha(vcf_path.read_bytes())
            new_format_path=roster_path.parent/'source-formats'/(case_id+'-small.json');_write(new_format_path,_bytes(source_format))
            case['source_format_bindings']=[new_format_path.relative_to(repository).as_posix()]
            case['annotation_expectations']=expectations
            case['regression_context']={'profile':suffix,'identifiers':'source-native-public-context','mechanisms':'source-required; unavailable remains unknown','observations':'plain-positive-HPO; no inferred negatives'}
            roster['cases'].append(case);truth['cases'].append({'truth_id':case['truth_id'],'case_id':case_id,'scope':case['scope'],'target_kind':'full_pair','members':members,'audience':'evaluator_only'})
            annotations['expectations'].extend(expectations)
    annotations['limitations'].append('New model-state interval expectations require actual VEP verification before empirical release acceptance; no synthetic annotation output is substituted.')
    _write(truth_path,_bytes(truth));_write(annotations_path,_bytes(annotations))
    roster['causal_truth'].update(_binding(repository,truth_path));roster['annotation_expectations'].update(_binding(repository,annotations_path))
    roster['limitations'].extend([
        'Original 22 fixtures and generator identities are retained; their previous held-out results are regression reuse, not renewed independence evidence.',
        'Fresh held-out model-state probes use generator seed 303 and a different public gene/chromosome from development; all related allele/phase perturbations share one split.',
        'Synthetic pair truth denotes a full-genotype recovery stress target, including the contradicted cis probe; it is not a clinical assertion or expected eligible submission.',
        'Current package extends original roster '+original_roster_sha+' without rebinding any historical result.',
    ])
    _write(roster_path,_bytes(roster))
    policy_path=repository/'config/track1/public-ranking-development-policy-family.json';policy=json.loads(policy_path.read_bytes());fixed=policy['fixed_settings']
    fixed['clinical_confidences']=list(REVIEW_STATUSES)
    fixed['population']={
        'derivation_sha256':_sha(derivation_bytes),
        'frequency_rule':{'rule_id':'development-selected-genomes-adjusted-AF','subset':'genomes','stratum':'adjusted','ancestry':None,'metric':'AF','threshold':derivation['selected_frequency_threshold'],'minimum_an':derivation['minimum_observed_development_an'],'allowed_filters':['PASS'],'region_flags':[],'unreliable_signals':['counts_missing','filters_unknown','inconsistent_ac_an','inconsistent_genotype_counts','low_an','source_filtered','zero_an']},
        'transform':{'kind':'capped_linear_excess','scale':'1','cap':'1'},
        'rescues':{'pathogenic_significances':fixed['clinical_positive_significances'],'pathogenic_review_statuses':fixed['clinical_review_statuses'],'pathogenic_conflicts':'block_assertion_rescue','pathogenic_match':'exact_native_condition_and_gene_target','technical_signals':['not_callable'],'technical_scope':'every_reliable_excess_record_or_every_unreliable_excess_record','founder':'exact_public_curated_exception','reduced_penetrance':'exact_public_curated_exception','combination':'any_applies'},
    }
    policy['resource_query']['population']='actual pinned joint BCF/CSI and both all-sites AN products; fixed development-derived rule instantiated for authenticated native contexts'
    policy['resource_query']['clinical']='actual Open Targets variant/evidence_eva and ClinVar VCF aggregates; only explicit native condition/target links may qualify'
    policy['limitations'].append('The added population separator is development-source-derived for synthetic stress assessment; it is not clinical maximum credible AF calibration.')
    _write(policy_path,_bytes(policy))
    _validate_materialized_call_inputs(
        repository, reference, contract, source_aliases=source_aliases,
    )


def _verification_cases(repository, roster_path, case_ids):
    from genosapiens.bounded_json import load_bounded_json
    roster, _content = load_bounded_json(repository / roster_path,
        reference_root=repository,maximum_bytes=2 * 1024 * 1024)
    if case_ids is None:
        selected = tuple(case for case in roster['cases'] if case['generator']['version'] == '2')
    else:
        names = tuple(case_ids)
        if not names or len(names) != len(set(names)):
            raise ValueError('annotation verification selection is invalid')
        index = {case['case_id']: case for case in roster['cases']}
        if len(index) != len(roster['cases']) or not set(names) <= index.keys():
            raise ValueError('annotation verification selection is invalid')
        selected = tuple(index[name] for name in names)
    if Path(roster_path).name == 'submission-population-source-roster.json':
        from genosapiens.track1_roster import POPULATION_GENERATOR_CASES
        if (case_ids is None or {case['case_id'] for case in selected} != set(POPULATION_GENERATOR_CASES)
                or any(case['generator']['version'] != '4' for case in selected)):
            raise ValueError('population annotation verification selection is invalid')
    return roster, selected


def _validate_materialized_call_inputs(repository, reference, contract, *,
                                      source_aliases, roster_path=ROSTER_PATH,
                                      case_ids=None):
    """Fail fixture publication on parser errors, before costly source replay.

    This is a provisional schema/REF probe, not full resource authentication or
    annotation evidence. The actual VEP verifier freshly authenticates the bundle.
    """
    from genosapiens.track1_public_validation import _prepare_call_input
    from genosapiens.admission_inventory import build_admission_inventory
    resources = SimpleNamespace(reference=reference, reference_contract=contract,
                                source_to_reference_contigs=source_aliases)
    _roster, selected = _verification_cases(repository, roster_path, case_ids)
    for case in selected:
        path = repository / next(value['path'] for value in case['runner_inputs'] if value['role'] == 'small_variant_vcf')
        inputs, _diagnostic = _prepare_call_input(path, (repository / case['source_format_bindings'][0]).read_bytes(), resources)
        build_admission_inventory(inputs)


def _annotation_verification_resources(repository, references, scratch, *,
                                       roster_path=None, case_ids=None,
                                       reuse_published_validation=False):
    """Authenticate only the resources used by this annotation-only probe.

    The full ranking runner independently authenticates phenotype, disease,
    clinical and population sources. This materializer neither ranks nor consumes
    those stores when replacing interval expectations from actual VEP output.
    """
    from genosapiens.reference_concordance import (
        load_construction_report, open_installed_reference,
    )
    from genosapiens.track1_public_validation import (
        REFERENCE_IDENTITY_PATH, REFERENCE_MANIFEST_SHA256,
        OPEN_TARGETS_CONTRACT_PATH, OPEN_TARGETS_MANIFEST_SHA256,
        VEP_CONTRACT_PATH, VEP_BUNDLE_SHA256,
        _authenticated_public_reference_aliases, _bound_bytes,
    )
    from genosapiens.vep_annotation import (
        OpenTargetsConstraintReference, VepReferenceCompatibility,
    )
    if any(name != 'lo' for _index, name in socket.if_nameindex()):
        raise ValueError
    if not scratch.is_relative_to(references / 'scratch/public-build'):
        raise ValueError
    if roster_path is not None:
        _verification_cases(repository, roster_path, case_ids)
    elif case_ids is not None:
        raise ValueError('annotation verification package is missing')
    identity_path = repository / REFERENCE_IDENTITY_PATH
    identity_sha256 = _sha(identity_path.read_bytes())
    contract = load_reference_identity(identity_path)
    manifest_path = references / 'manifests/grch38/challenge-grch38-masked-v2-plus-hs38d1.a92c50f61c1b.json'
    _bound_bytes(manifest_path, REFERENCE_MANIFEST_SHA256)
    with open_installed_reference(
        contract, identity_path=identity_path,
        reference_manifest_path=manifest_path, reference_root=references,
        reuse_published_validation=reuse_published_validation,
    ) as (manifest_sha256, _fetch, contigs):
        if manifest_sha256 != REFERENCE_MANIFEST_SHA256 or len(contigs) != contract.sequence_count:
            raise ValueError
    role = contract.reference_role('public_annotation_benchmark')
    reference = IndexedFasta(
        references / contract.expected_output(role.fasta_output_role).path,
        references / contract.expected_output(role.fai_output_role).path,
    )
    construction = load_construction_report(references / contract.expected_output('construction_report').path)
    if (not construction.accepted or construction.public_annotation_benchmark_identity_id != role.identity_id
            or construction.total_records != len(reference.index)):
        raise ValueError
    ot_manifest = references / 'manifests/open-targets/open-targets-platform-26.06.0cae4540cfea.json'
    _bound_bytes(ot_manifest, OPEN_TARGETS_MANIFEST_SHA256)
    _bound_bytes(references / 'manifests/ensembl/vep/ensembl-vep-116-grch38.f1c47e44cd07.bundle.json', VEP_BUNDLE_SHA256)
    if not (repository / VEP_CONTRACT_PATH).is_file():
        raise ValueError
    source_aliases = _authenticated_public_reference_aliases(
        identity_path=identity_path,
        reference_manifest_path=manifest_path,
        reference_root=references,
        reference_contract=contract,
        reference=reference,
    )
    return SimpleNamespace(
        repository_root=repository, reference_root=references, scratch_root=scratch,
        reference_contract=contract, reference=reference,
        source_to_reference_contigs=source_aliases,
        vep_constraint=OpenTargetsConstraintReference(
            repository / OPEN_TARGETS_CONTRACT_PATH,
            references / 'releases/open-targets/26.06', ot_manifest,
            OPEN_TARGETS_MANIFEST_SHA256,
        ),
        vep_compatibility=VepReferenceCompatibility(
            identity_path, identity_sha256, references, manifest_path,
            REFERENCE_MANIFEST_SHA256,
        ),
    )


def verify_model_annotations(repository, references, scratch, *,
                             roster_path=ROSTER_PATH, case_ids=None):
    """Replace planned interval expectations only with actual VEP readback."""
    from genosapiens.admission_inventory import build_admission_inventory
    from genosapiens.annotation_admission import replay_annotation_admission
    from genosapiens.evidence_model import AlleleEvidence
    from genosapiens.evidence_serialization import from_json
    from genosapiens.track1_public_validation import (
        _annotation_artifacts, _prepare_call_input,
    )
    from genosapiens.track1_validation_runner import _member_keys
    from genosapiens.vep_annotation import VepRequest
    print('model-state annotation verification: reference authentication', flush=True)
    resources=_annotation_verification_resources(repository,references,scratch,
        roster_path=roster_path,case_ids=case_ids)
    print('model-state annotation verification: request preparation', flush=True)
    roster,selected=_verification_cases(repository,roster_path,case_ids)
    roster_path=repository/roster_path
    annotations_path=repository/roster['annotation_expectations']['path'];annotations=json.loads(annotations_path.read_bytes())
    requests=[]
    for case in selected:
        vcf_path=repository/next(item['path'] for item in case['runner_inputs'] if item['role']=='small_variant_vcf')
        source_format=(repository/case['source_format_bindings'][0]).read_bytes()
        inputs,_diagnostic=_prepare_call_input(vcf_path,source_format,resources)
        # Sequence-derived allele IDs recur across perturbations. Keep every
        # source call and its provenance, giving only the enclosing request a
        # source-input namespace, as the source-set inventory does. Native and
        # submission identities, sequence, genotype and annotation are unchanged.
        requests.extend(VepRequest(replace(
            allele, allele_id='annotation-probe-' + _sha(
                (inputs.input_sha256 + '\n' + allele.allele_id).encode('ascii')
            )[:32],
        )) for allele in build_admission_inventory(inputs).alleles)
    requests=tuple(sorted(requests,key=lambda value:value.allele.allele_id))
    print('model-state annotation verification: actual VEP execution', flush=True)
    namespace='model-roster-annotation-verification-' + _sha(_bytes([
        _sha(roster_path.read_bytes()), [case['case_id'] for case in selected]]))[:16]
    artifacts=_annotation_artifacts(requests,resources=resources,publication_id=namespace)
    print('model-state annotation verification: annotation admission replay', flush=True)
    replay_annotation_admission(artifacts,expected_binding_sha256=_sha(artifacts['binding.json']),requests=requests)
    actual={}
    for value in json.loads(artifacts['alleles.json']):
        allele=from_json(json.dumps(value),AlleleEvidence)
        for member in _member_keys(allele):
            genes=sorted({t.hgnc_symbol for t in allele.transcript_consequences or () if t.hgnc_symbol})
            consequences=sorted({term for t in allele.transcript_consequences or () for term in t.consequences or () if term})
            entry={'member':member,'gene_symbols':genes,'most_severe_consequence':None,'consequences':consequences,'assessment_status':'actual-VEP-116-readback-and-admission-replay; no invented severity selection'}
            if member in actual and actual[member]!=entry:raise ValueError
            actual[member]=entry
    refreshed=[]
    selected_ids={case['case_id'] for case in selected}
    for case in roster['cases']:
        if case['case_id'] in selected_ids:
            case['annotation_expectations']=[{**actual[entry['member']],'split':case['split']} for entry in case['annotation_expectations']]
            case['metric_applicability']['gene_rank']=any(entry['gene_symbols'] for entry in case['annotation_expectations'])
        refreshed.extend(case['annotation_expectations'])
    annotations['expectations']=refreshed
    print('model-state annotation verification: expectation publication', flush=True)
    annotations['source_runs']['current-model-state']={'kind':'current-model-state-actual-VEP', 'binding_sha256':_sha(artifacts['binding.json']),'request_count':len(requests),'artifact_namespace':namespace}
    annotations['limitations'].append('New model-state expectations were replaced from actual VEP readback; source interval selection alone did not establish gene overlap or clinical qualification.')
    _write(annotations_path,_bytes(annotations));roster['annotation_expectations'].update(_binding(repository,annotations_path));_write(roster_path,_bytes(roster))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repository-root',type=Path,required=True);parser.add_argument('--reference-root',type=Path,required=True)
    parser.add_argument('--verify-annotations',action='store_true');parser.add_argument('--scratch-root',type=Path)
    parser.add_argument('--roster-path',type=Path,default=ROSTER_PATH)
    parser.add_argument('--case-id',action='append',dest='case_ids')
    args=parser.parse_args()
    try:
        repository=args.repository_root.resolve(strict=True);references=args.reference_root.resolve(strict=True)
        if args.verify_annotations:
            if args.scratch_root is None:raise ValueError
            verify_model_annotations(repository,references,args.scratch_root.resolve(strict=True),
                roster_path=args.roster_path,case_ids=args.case_ids)
            print('public model-state annotations verified with actual VEP and replay')
        else:
            materialize_model_state_roster(repository,references)
            print('public model-state roster materialized; actual-tool verification remains required')
    except Exception:
        print('public model-state materialization failed; no scientific completion claimed')
        raise SystemExit(1) from None


if __name__=='__main__':main()
