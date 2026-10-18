#!/usr/bin/env python3
"""Check the publication without MATLAB or third-party Python packages.

Checks tracked files in a Git clone, or the exported tree otherwise. Generated
results are not scientific reference inputs and must remain outside Git.
"""
import argparse
import ast
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
REQUIRED = {
    "README.md", "LICENSE", "CONTRIBUTING.md", "CHANGELOG.md", "CITATION.cff",
    ".github/workflows/verify.yml", "data/manifest.json", "docs/index.md",
    "docs/getting-started.md", "docs/architecture.md", "docs/solid-fem.md",
    "docs/api-reference.md", "docs/examples.md", "docs/data-and-preprocessing.md",
    "docs/verification.md", "docs/reduced-model.md", "docs/experimental-data.md",
    "docs/references.md", "docs/roadmap.md", "matlab/run_verification.m",
    "matlab/run_fem_demo.m", "matlab/run_fem_reference.m", "matlab/run_fem_convergence.m",
    "matlab/tests/run_tests.m", "matlab/tests/run_fem_tests.m",
    "matlab/tests/run_tail_tests.m", "matlab/run_tail_benchmarks.m", "matlab/run_tail_response.m",
    "docs/trailing-edge.md", "docs/case-catalog.md", "THIRD_PARTY_NOTICES.md",
}


def published_paths(include_working_tree=False):
    try:
        top = subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "--show-toplevel"],
            stderr=subprocess.DEVNULL, text=True).strip()
        if Path(top).resolve() == ROOT:
            raw = subprocess.check_output(["git", "-C", str(ROOT), "ls-files", "-z"] +
                                          (["--cached", "--others", "--exclude-standard"] if include_working_tree else []))
            paths = [ROOT / name for name in raw.decode().split("\0") if name]
            if paths:
                return paths
    except (OSError, subprocess.CalledProcessError):
        pass
    ignored = {".git", ".cache", ".venv", "results", "__pycache__", ".publication-check"}
    return [p for p in ROOT.rglob("*") if p.is_file()
            and not any(part in ignored for part in p.relative_to(ROOT).parts)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-working-tree", action="store_true",
                        help="Also validate untracked, non-ignored new deliverables without staging them")
    args = parser.parse_args()
    errors = []
    paths = published_paths(args.include_working_tree)
    names = {str(p.relative_to(ROOT).as_posix()) for p in paths}
    for name in sorted(REQUIRED - names):
        errors.append(f"Missing required publication file: {name}")

    forbidden_parts = {"results", ".cache", ".venv", "__pycache__"}
    forbidden_suffixes = {".mat", ".fig", ".zip", ".stp", ".step", ".ipt", ".stl", ".pyc"}
    for p in paths:
        name = p.relative_to(ROOT).as_posix()
        if any(part in forbidden_parts for part in p.relative_to(ROOT).parts):
            errors.append(f"Generated/environment content tracked: {name}")
        if p.suffix.lower() in forbidden_suffixes or name.startswith("reference_preview_"):
            errors.append(f"Non-runtime original/export tracked: {name}")
        if not p.is_file() or p.is_symlink():
            errors.append(f"Missing or symlinked publication file: {name}")
            continue
        if p.stat().st_size > 20 * 1024 * 1024:
            errors.append(f"Unexpected individual file size: {name}")
        if p.suffix == ".json":
            try:
                json.loads(p.read_text())
            except (ValueError, OSError) as exc:
                errors.append(f"Invalid JSON: {name}: {exc}")
        if p.suffix == ".py":
            try:
                ast.parse(p.read_text(), filename=name)
            except SyntaxError as exc:
                errors.append(f"Invalid Python: {name}: {exc}")

    manifest = json.loads((ROOT / "data/manifest.json").read_text())
    for entry in manifest["files"]:
        p = ROOT / entry["path"]
        if entry["path"] not in names:
            errors.append(f"Manifest input is absent from publication: {entry['path']}")
            continue
        if hashlib.sha256(p.read_bytes()).hexdigest() != entry["sha256"]:
            errors.append(f"Input checksum mismatch: {entry['path']}")
        if p.suffix == ".csv":
            with p.open(newline="") as handle:
                rows = csv.reader(handle)
                header = next(rows)
                count = 0
                for row in rows:
                    count += 1
                    if len(row) != len(header):
                        errors.append(f"Nonrectangular CSV: {entry['path']} row {count+1}")
                        break
            if count != entry["rows"] or len(header) != entry["columns"]:
                errors.append(f"Input table dimensions differ: {entry['path']}")
    data_names = {name for name in names if name.startswith("data/")}
    expected_data = {entry["path"] for entry in manifest["files"]} | {"data/README.md", "data/manifest.json"}
    for name in sorted(data_names - expected_data):
        errors.append(f"Unexpected data file not in runtime manifest: {name}")

    source_hash = manifest["source_step_sha256"]
    for directory in ("reference", "fem_reference", "fem_refined"):
        meta_name = "cad_info.json" if directory == "reference" else "mesh_info.json"
        meta = json.loads((ROOT / "data" / directory / meta_name).read_text())
        if meta["units"] != "m" or meta["source_sha256"] != source_hash:
            errors.append(f"Inconsistent units/source: {directory}")
    audit = json.loads((ROOT / "data/reference/cad_volume_audit.json").read_text())
    if not audit["brep_valid"] or audit["source_sha256"] != source_hash:
        errors.append("Invalid or mismatched stored BREP audit")

    links = 0
    for p in paths:
        if p.suffix != ".md":
            continue
        text = p.read_text()
        text = re.sub(r"```[\s\S]*?```", "", text)
        if text.count("```") % 2:
            errors.append(f"Unclosed code fence: {p.relative_to(ROOT)}")
        for match in re.finditer(r"\[[^\]\n]*\]\(([^\s)]+)(?:\s+\"[^\"]*\")?\)", text):
            target = match.group(1).strip("<>")
            url = urlsplit(target)
            if url.scheme or url.netloc or not url.path:
                continue
            links += 1
            resolved = (p.parent / unquote(url.path)).resolve()
            try:
                relative = resolved.relative_to(ROOT).as_posix()
            except ValueError:
                errors.append(f"Documentation link escapes repository: {target}")
                continue
            if relative not in names and not resolved.is_dir():
                errors.append(f"Broken/unpublished local link in {p.relative_to(ROOT)}: {target}")

    if errors:
        for error in errors:
            print("ERROR:", error, file=sys.stderr)
        raise SystemExit(1)
    print(f"PASS: {len(names)} publication files, {len(manifest['files'])} hashed inputs, {links} local documentation links.")
    print("Numerical verification is separate: addpath('matlab'); run_verification();")


if __name__ == "__main__":
    main()
