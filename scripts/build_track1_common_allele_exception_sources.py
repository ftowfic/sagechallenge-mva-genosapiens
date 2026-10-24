#!/usr/bin/env python3
"""Freeze public ClinGen exception inputs without case or truth data."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from urllib.request import Request, urlopen


BA1_URL = "https://clinicalgenome.org/docs/ba1-exception-list/"
BA1_RESOLVED_URL = "https://clinicalgenome.org/site/assets/files/3460/ba1_exception_list_07_30_2018.pdf"
REGISTRY_URL = "https://cspec.genome.network/cspec/api/svis"
TERMS_URL = "https://www.clinicalgenome.org/docs/terms-of-use/"
DOCUMENT_API = "https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/{svi_id}"
POLICY_PATH = "config/track1/public-ranking-development-policy-family.json"
MAX_SOURCE_BYTES = 64 * 1024 * 1024
ID = re.compile(r"^GN[0-9]{3}$")


class SourceSnapshotError(RuntimeError):
    pass


def fail() -> None:
    raise SourceSnapshotError("public common-allele exception source snapshot failed")


def canonical(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def fetch(url: str) -> tuple[str, bytes]:
    request = Request(url, headers={"User-Agent": "genosapiens-public-resource-freeze/1"})
    with urlopen(request, timeout=120) as response:
        content = response.read(MAX_SOURCE_BYTES + 1)
        final_url = response.geturl()
    if not content or len(content) > MAX_SOURCE_BYTES:
        fail()
    return final_url, content


def fetch_or_existing(path: Path, url: str) -> tuple[str, bytes]:
    """Resume an interrupted immutable snapshot without refetching its prefix."""

    if path.exists():
        if path.is_symlink() or not path.is_file():
            fail()
        content = path.read_bytes()
        if not content or len(content) > MAX_SOURCE_BYTES:
            fail()
        return (BA1_RESOLVED_URL if url == BA1_URL else url), content
    return fetch(url)


def atomic(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".partial-", dir=path.parent)
    try:
        os.fchmod(descriptor, 0o644)
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    except FileExistsError:
        if path.read_bytes() != content:
            fail()
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def artifact(role: str, requested_url: str, final_url: str, relative: str, content: bytes) -> dict:
    return {
        "role": role,
        "requested_url": requested_url,
        "resolved_url": final_url,
        "path": relative,
        "size_bytes": len(content),
        "sha256": digest(content),
    }


def checked_json(content: bytes) -> dict:
    try:
        value = json.loads(content)
    except (ValueError, TypeError, UnicodeError):
        fail()
    if type(value) is not dict:
        fail()
    return value


def identifier(entry: dict) -> str:
    value = entry.get("@id")
    if type(value) is not str:
        fail()
    svi_id = value.rsplit("/", 1)[-1]
    if ID.fullmatch(svi_id) is None:
        fail()
    return svi_id


def gene_rows(rule_sets: list) -> list[dict]:
    genes: dict[tuple[str, str], dict] = {}
    for rule_set in rule_sets:
        if type(rule_set) is not dict:
            fail()
        listed = rule_set.get("genes")
        if listed is None:
            continue
        if type(listed) is not list:
            fail()
        for gene in listed:
            if type(gene) is not dict or type(gene.get("label")) is not str:
                fail()
            source = gene.get("@id")
            if type(source) is not str:
                fail()
            match = re.search(r"HGNC:[0-9]+", source)
            hgnc_id = match.group(0) if match else None
            key = (gene["label"], source)
            genes[key] = {
                "symbol": gene["label"],
                "hgnc_id": hgnc_id,
                "source_id": source,
            }
    return [genes[key] for key in sorted(genes)]


def ba1_rows(rule_sets: list) -> list[dict]:
    rows = []
    for rule_set in rule_sets:
        rule_set_id = rule_set.get("@id")
        criteria = rule_set.get("criteriaCodes")
        if type(rule_set_id) is not str or type(criteria) is not list:
            fail()
        matches = [item for item in criteria if type(item) is dict and item.get("label") == "BA1"]
        if len(matches) > 1:
            fail()
        for criterion in matches:
            strengths = criterion.get("evidenceStrengths")
            if type(strengths) is not list:
                fail()
            standalone = [item for item in strengths if type(item) is dict
                          and item.get("label") == "Stand Alone"]
            if len(standalone) != 1:
                fail()
            strength = standalone[0]
            rows.append({
                "rule_set_id": rule_set_id,
                "criteria_id": criterion.get("@id"),
                "applicability": strength.get("applicability"),
                "description": strength.get("description"),
            })
    return rows


def build(reference_root: Path, repository_root: Path, retrieved_date: str) -> Path:
    if re.fullmatch(r"20[0-9]{2}-[0-9]{2}-[0-9]{2}", retrieved_date) is None:
        fail()
    release_relative = f"releases/track1-common-allele-exceptions/{retrieved_date}/sources"
    release = reference_root / release_relative
    manifest_dir = reference_root / "manifests/track1-common-allele-exceptions"
    if release.is_symlink() or manifest_dir.is_symlink():
        fail()

    sources = []
    for role, url, filename in (
        ("clingen_ba1_exception_list", BA1_URL, "clingen-ba1-exception-list-2018-07-30.pdf"),
        ("clingen_terms", TERMS_URL, "clingen-terms.html"),
        ("cspec_registry_listing", REGISTRY_URL, "cspec-registry-listing.json"),
    ):
        relative = f"{release_relative}/{filename}"
        final_url, content = fetch_or_existing(reference_root / relative, url)
        atomic(reference_root / relative, content)
        sources.append(artifact(role, url, final_url, relative, content))

    listing_content = (reference_root / sources[-1]["path"]).read_bytes()
    listing = checked_json(listing_content)
    entries = listing.get("data")
    if type(entries) is not list or not entries:
        fail()
    ids = [identifier(entry) for entry in entries]
    if len(ids) != len(set(ids)):
        fail()
    released = sorted((entry for entry in entries if entry.get("status") == "Released"),
                      key=identifier)
    documents = []
    for entry in released:
        svi_id = identifier(entry)
        requested_url = DOCUMENT_API.format(svi_id=svi_id)
        relative = f"{release_relative}/cspec/{svi_id}.json"
        final_url, content = fetch_or_existing(reference_root / relative, requested_url)
        document = checked_json(content)
        rule_sets = document.get("ruleSets")
        if (type(rule_sets) is not list or not rule_sets
                or identifier(document) != svi_id
                or document.get("cspecStatus") != "Released"
                or document.get("version") != entry.get("version")):
            fail()
        atomic(reference_root / relative, content)
        ba1 = ba1_rows(rule_sets)
        genes = gene_rows(entry.get("ruleSets"))
        documents.append({
            "svi_id": svi_id,
            "version": document["version"],
            "last_updated": document.get("lastUpdated"),
            "ui_url": entry.get("url"),
            "api_url": final_url,
            "path": relative,
            "size_bytes": len(content),
            "sha256": digest(content),
            "rule_set_count": len(rule_sets),
            "genes": genes,
            "gene_coverage": "present" if genes else "missing_declared",
            "ba1_coverage": "present" if ba1 else "missing_declared",
            "ba1_records": ba1,
            "ba1_records_sha256": digest(canonical(ba1)),
        })

    if {item["svi_id"] for item in documents} != {identifier(item) for item in released}:
        fail()
    policy = repository_root / POLICY_PATH
    policy_content = policy.read_bytes()
    downstream = {
        "role": "selected_downstream_public_exception_curation",
        "policy_path": POLICY_PATH,
        "policy_sha256": digest(policy_content),
        "founder_records": [],
        "reduced_penetrance_records": [],
        "hypomorphic_records": [],
        "empty_source_declaration": True,
        "reason": "selected public downstream policy defines rescue adapters but freezes no curated exception records",
    }
    manifest = {
        "schema": "genosapiens.track1-common-allele-exception-source-snapshot",
        "version": 1,
        "classification": "public",
        "snapshot_id": f"track1-common-allele-exception-sources-{retrieved_date}",
        "retrieved_date": retrieved_date,
        "license": "CC0-1.0",
        "license_url": "https://creativecommons.org/publicdomain/zero/1.0/",
        "terms_url": TERMS_URL,
        "challenge_or_truth_inputs_used": False,
        "registry_listing_count": len(entries),
        "registry_released_count": len(released),
        "sources": sources,
        "registry_documents": documents,
        "downstream_public_curation": downstream,
    }
    content = canonical(manifest)
    path = manifest_dir / f"source-snapshot.{digest(content)[:12]}.json"
    atomic(path, content)
    if path.read_bytes() != content:
        fail()
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--retrieved-date", required=True)
    arguments = parser.parse_args()
    path = build(arguments.reference_root.resolve(), arguments.repository_root.resolve(),
                 arguments.retrieved_date)
    print(path)


if __name__ == "__main__":
    main()
