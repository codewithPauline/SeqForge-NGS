#!/usr/bin/env python3
"""Verify a prepared benchmark bundle before running the caller."""
import argparse
import json
from pathlib import Path
from benchmark_io import verify

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bundle')
    args = parser.parse_args()
    try:
        Path('verification.json').write_text(json.dumps(verify(args.bundle), indent=2) + '\n')
    except (OSError, KeyError, ValueError) as exc:
        parser.exit(2, f'Benchmark verification failed: {exc}\n')
