#!/usr/bin/env python3
"""Public-only acquisition and offline discovery catalog preparation."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT))


def ensure_discovery_resources(repository,reference_root,*,memory_gib=4):
    from genosapiens.track2_discovery_resources import (acquire_sources,load_discovery_resources,
        public_path,file_digest,file_state,publish_catalog)
    from genosapiens.track2_run_state import atomic_record
    from genosapiens.track2_published_genotype import canonical,sha256
    from genosapiens.track2_construction_dependencies import python_dependencies
    from scripts.run_track2 import Docker,mount
    repository,reference_root=Path(repository).absolute(),Path(reference_root).absolute()
    read=lambda path:json.loads(Path(path).read_bytes())
    contract=read(repository/'config/track2/discovery-sources.json')
    def tissue_record(record):
        from genosapiens.track2_tissue_resources import prepare_tissue_companion
        descriptor=read(public_path(reference_root,record['path']))
        if file_digest(public_path(reference_root,record['path']))!=record['sha256']:
            raise ValueError('The retained discovery descriptor changed.')
        if descriptor['schema']=='genosapiens.track2-discovery-resources/v2':
            load_discovery_resources(descriptor,reference_root).close()
            return record
        tissue=prepare_tissue_companion(descriptor,reference_root,contract,
            read(repository/'config/track2/exploratory-policy-v1.json'),
            read(repository/'config/track2/exploratory-policy-v2.json'))
        path=public_path(reference_root,tissue['tissue_companion']['path']).parent/'manifest.json'
        installed=dict(schema='genosapiens.track2-discovery-installation/v2',
            path=path.relative_to(reference_root).as_posix(),sha256=file_digest(path),build_id=tissue['build_id'])
        atomic_record(reference_root/'derived/track2/discovery/current-tissue.json',installed)
        return installed
    # Query-only preparation uses retained metadata before any builder request.
    # Keep the original current.json and its producer completion unchanged.
    retained=reference_root/'derived/track2/discovery/current.json'
    if retained.exists():
        return tissue_record(read(retained))
    from genosapiens.track2_discovery_resources import require_discovery_admission
    input_admission=require_discovery_admission(contract, reference_root)
    drug_contract=read(repository/contract['open_targets']['contract_path'])
    from genosapiens.portable_resources import prepare_open_targets_metadata
    prepare_open_targets_metadata(reference_root, drug_contract)
    builder=read(repository/'config/track2/discovery-builder-runtime.json')
    subprocess.run(['docker','image','inspect',builder['image']],stdout=subprocess.DEVNULL,check=True)
    root=public_path(reference_root,'derived/track2/discovery');root.mkdir(parents=True,exist_ok=True)
    with (root/'.prepare.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        acquired=acquire_sources(drug_contract,reference_root)
        names=[contract['open_targets']['target_manifest_path']]
        for source,key in [('lincs','metadata_manifest'),('lincs','matrix_receipt'),('lincs','transfer_receipt'),('fda','manifest')]:
            names.append(contract[source]['root']+'/'+contract[source][key])
        states={}
        for source in ['fda','lincs']:
            selected=contract[source]['files' if source=='fda' else 'metadata_files']
            for filename in selected.values():
                path=public_path(reference_root,contract[source]['root']+'/'+filename)
                if not path.is_file():raise ValueError('Restore the declared public source asset: '+str(path))
                states[path.relative_to(reference_root).as_posix()]=file_state(path)
        code_names=['track2_discovery_resources.py','track2_discovery_catalog.py','track2_medication_catalog.py','track2_lincs_selection.py']
        code_files={name:file_digest(repository/'src/genosapiens'/name) for name in code_names}
        request=sha256(canonical(dict(contract=contract,drug_contract=drug_contract,acquisition_sha256=acquired['sha256'],code_files=code_files,
            input_admission=input_admission,
            policy=read(repository/'config/track2/exploratory-policy-v1.json'),builder=builder,
            manifests={name:file_digest(public_path(reference_root,name)) for name in names},metadata_states=states,
            code=python_dependencies(repository,[('genosapiens.track2_discovery_resources','build_catalog')]))))
        current=root/'current.json'
        if current.exists():
            record=read(current)
            if record['request_sha256']==request:
                descriptor=public_path(reference_root,record['path'])
                if file_digest(descriptor)!=record['sha256']:raise ValueError('The current discovery descriptor changed.')
                load_discovery_resources(descriptor,reference_root).close()
                return tissue_record(record)
        staging=root/('.building-'+request)
        if staging.exists():
            from genosapiens.track2_run_state import files
            files(staging);shutil.rmtree(staging)
        staging.mkdir(mode=0o700)
        atomic_record(root/'preparation-progress.json',dict(status='running',request_sha256=request))
        docker=Docker(builder['image'],memory_gib)
        snapshot=tempfile.TemporaryDirectory(prefix='.builder-source-',dir=root)
        frozen=Path(snapshot.name)
        for name in ['src','config','schemas']:
            shutil.copytree(repository/name,frozen/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        if code_files!={name:file_digest(frozen/'src/genosapiens'/name) for name in code_names}:
            snapshot.cleanup()
            raise ValueError('Discovery builder sources changed during preparation. Repeat public preparation.')
        options=['--env','PYTHONPATH='+str(repository/'src'),'--workdir',str(repository)]
        for name in ['src','config','schemas']:
            options+=['--mount',mount(frozen/name,repository/name)]
        for name in ['releases','manifests','lincs','fda']:
            options+=['--mount',mount(reference_root/name)]
        # Admission can reuse an existing public catalog when no local seal exists.
        # The worker needs the same catalog evidence that the host accepted.
        prior_catalogs=public_path(reference_root,'derived/track2/discovery/v1')
        if prior_catalogs.is_dir():
            options+=['--mount',mount(prior_catalogs)]
        options+=['--mount',mount(staging,readonly=False)]
        code='''import json,sys
from pathlib import Path
from genosapiens.track2_discovery_resources import build_catalog
build_catalog(json.loads(Path(sys.argv[1]).read_bytes()),sys.argv[2],sys.argv[3])
'''
        try:
            status=docker.run(options,['-c',code,str(repository/'config/track2/discovery-sources.json'),str(reference_root),str(staging)],
                log=root/'builder.log')
            if status:raise RuntimeError('The public discovery builder failed. Its outputs remain incomplete.')
            manifest=read(staging/'manifest.json')
            target=root/'v1'/manifest['build_id']
            if target.exists():
                load_discovery_resources(target/'manifest.json',reference_root).close()
                shutil.rmtree(staging)
                descriptor=target/'manifest.json'
            else:descriptor=publish_catalog(staging,target)
            record=dict(schema='genosapiens.track2-discovery-installation/v1',request_sha256=request,
                path=descriptor.relative_to(reference_root).as_posix(),sha256=file_digest(descriptor),build_id=manifest['build_id'])
            atomic_record(current,record)
            atomic_record(root/'preparation-progress.json',dict(status='complete',**record))
            return tissue_record(record)
        except BaseException:
            docker.cleanup();atomic_record(root/'preparation-progress.json',dict(status='incomplete',request_sha256=request));raise
        finally:snapshot.cleanup()


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-root',type=Path,default=Path('/mnt/reference'))
    parser.add_argument('--memory-gib',type=int,default=4)
    args=parser.parse_args(argv)
    print(json.dumps(ensure_discovery_resources(ROOT,args.reference_root,memory_gib=args.memory_gib)))


if __name__=='__main__':main()
