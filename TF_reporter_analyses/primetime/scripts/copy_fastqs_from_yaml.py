#!/usr/bin/env python3
"""
Scan all .yaml files in a primetime folder, extract FASTQ paths, and copy them to a destination folder.

Features:
- Tries to use PyYAML to parse YAML files. If PyYAML is not available, falls back to a regex-based extraction for FASTQ paths.
- Supports dry-run, verbose, preserving directory structure under destination, and skipping existing files.

Usage examples:
  python copy_fastqs_from_yaml.py --primetime /path/to/primetime --dest /path/to/collected_fastqs --dry-run
  python copy_fastqs_from_yaml.py --dest ./collected_fastqs --preserve-structure

"""
import argparse
import os
import sys
import shutil
import glob
import re
from pathlib import Path


def find_yaml_files(primetime_dir):
    pattern = os.path.join(primetime_dir, "*.yaml")
    return sorted(glob.glob(pattern))


def extract_fastqs_with_pyyaml(path):
    try:
        import yaml
    except Exception:
        raise

    with open(path, "r") as fh:
        data = yaml.safe_load(fh)

    found = []
    # YAML structure in these files places inputs under INPUT_DATA
    if not isinstance(data, dict):
        return found

    input_data = data.get("INPUT_DATA") or data.get("input_data")
    if not isinstance(input_data, dict):
        return found

    for key, val in input_data.items():
        if not isinstance(val, dict):
            continue
        fastq = val.get("fastq")
        if fastq is None:
            continue
        if isinstance(fastq, list):
            found.extend([str(x) for x in fastq if x])
        else:
            found.append(str(fastq))

    return found


def extract_fastqs_with_regex(path):
    # Fallback: find absolute or relative paths that end with .fastq or .fastq.gz
    fastq_re = re.compile(r"(/[^\s'\"]+?\.fastq(?:\.gz)?)")
    found = []
    with open(path, "r") as fh:
        for line in fh:
            for m in fastq_re.findall(line):
                found.append(m)
    return found


def collect_fastq_paths(yaml_files, verbose=False):
    paths = []
    pyyaml_available = True
    try:
        import yaml  # noqa: F401
    except Exception:
        pyyaml_available = False
        if verbose:
            print("PyYAML not available — falling back to regex extraction.")

    for y in yaml_files:
        if verbose:
            print(f"Parsing: {y}")
        try:
            if pyyaml_available:
                found = extract_fastqs_with_pyyaml(y)
            else:
                found = extract_fastqs_with_regex(y)
        except Exception as e:
            if verbose:
                print(f"  Failed parsing with PyYAML for {y}: {e}. Falling back to regex.")
            found = extract_fastqs_with_regex(y)

        for f in found:
            if f:
                paths.append(f)

    # normalize
    normed = []
    for p in paths:
        p = os.path.expanduser(p)
        normed.append(os.path.normpath(p))

    return sorted(set(normed))


def copy_files(paths, dest_dir, preserve_structure=False, dry_run=True, overwrite=False, verbose=False):
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    missing = []
    copied = []
    skipped = []

    for p in paths:
        src = Path(p)
        if not src.exists():
            missing.append(p)
            if verbose:
                print(f"MISSING: {p}")
            continue

        if preserve_structure:
            # create a path under dest that mirrors the absolute path without leading '/'
            rel = src.as_posix().lstrip("/")
            dst = dest_dir.joinpath(rel)
            dst.parent.mkdir(parents=True, exist_ok=True)
        else:
            dst = dest_dir.joinpath(src.name)

        if dst.exists() and not overwrite:
            skipped.append(str(dst))
            if verbose:
                print(f"SKIP (exists): {dst}")
            continue

        if dry_run:
            copied.append((str(src), str(dst)))
            if verbose:
                print(f"DRYRUN copy: {src} -> {dst}")
        else:
            shutil.copy2(src, dst)
            copied.append((str(src), str(dst)))
            if verbose:
                print(f"COPIED: {src} -> {dst}")

    return {
        "copied": copied,
        "skipped": skipped,
        "missing": missing,
    }


def main():
    parser = argparse.ArgumentParser(description="Collect FASTQ files referenced in primetime YAMLs and copy them to a destination folder.")
    parser.add_argument("--primetime", dest="primetime_dir", default=None,
                        help="Path to primetime folder (default: infer from script location)")
    parser.add_argument("--dest", required=True, help="Destination folder where FASTQ files will be copied")
    parser.add_argument("--preserve-structure", action="store_true", help="Preserve full path structure under destination (may create deep dirs)")
    parser.add_argument("--dry-run", action="store_true", default=False, help="Only show what would be copied")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files in destination")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    args = parser.parse_args()

    if args.primetime_dir:
        primetime_dir = args.primetime_dir
    else:
        # assume this script lives in primetime/scripts/
        script_dir = Path(__file__).resolve().parent
        primetime_dir = script_dir.parent.as_posix()

    if args.verbose:
        print(f"Primetime dir: {primetime_dir}")
        print(f"Destination: {args.dest}")

    yaml_files = find_yaml_files(primetime_dir)
    if not yaml_files:
        print(f"No YAML files found in {primetime_dir}")
        sys.exit(1)

    paths = collect_fastq_paths(yaml_files, verbose=args.verbose)
    if not paths:
        print("No FASTQ paths extracted from YAML files.")
        sys.exit(1)

    if args.verbose:
        print(f"Found {len(paths)} unique FASTQ paths")

    result = copy_files(paths, args.dest, preserve_structure=args.preserve_structure,
                        dry_run=args.dry_run, overwrite=args.overwrite, verbose=args.verbose)

    print("\nSummary:")
    print(f"  to-copy (count): {len(result['copied'])}")
    print(f"  skipped (exists): {len(result['skipped'])}")
    print(f"  missing: {len(result['missing'])}")

    if result['missing']:
        print("\nMissing files (not found):")
        for m in result['missing']:
            print("  ", m)

    if args.dry_run:
        print("\nDry-run mode: no files were actually copied. Re-run without --dry-run to perform copying.")


if __name__ == "__main__":
    main()
