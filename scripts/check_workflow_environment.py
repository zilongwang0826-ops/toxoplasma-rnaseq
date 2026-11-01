#!/usr/bin/env python3
"""Compare an exported conda explicit environment with its committed lock."""
import argparse
from pathlib import Path
import re
import sys
from urllib.parse import urlsplit


def explicit_packages(path):
    """Read a nonempty explicit list without discarding URL/hash differences."""
    packages = set()
    marker_seen = False
    for number, raw in enumerate(Path(path).read_text().splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        location = f'{path}:{number}'
        if line == '@EXPLICIT':
            if marker_seen or packages:
                raise ValueError(f'{location}: duplicate or misplaced @EXPLICIT marker')
            marker_seen = True
            continue
        if not marker_seen:
            raise ValueError(f'{location}: package entry before @EXPLICIT marker')
        try:
            url = urlsplit(line)
            url.port  # Validate an optional port without rewriting the URL.
        except ValueError as error:
            raise ValueError(f'{location}: malformed package URL') from error
        if (url.scheme != 'https' or not url.hostname or url.username is not None
                or url.password is not None or url.query
                or not url.path.endswith(('.conda', '.tar.bz2'))
                or any(character.isspace() for character in line)
                or ('#' in line and not re.fullmatch(r'(?:[0-9a-fA-F]{32}|[0-9a-fA-F]{64})', url.fragment))):
            raise ValueError(f'{location}: expected an HTTPS package URL with an optional MD5/SHA-256 hash')
        if line in packages:
            raise ValueError(f'{location}: duplicate package entry')
        packages.add(line)
    if not marker_seen or not packages:
        raise ValueError(f'{path}: explicit list must contain @EXPLICIT and at least one package')
    return packages


def compare_environment(lock, exported):
    expected = explicit_packages(lock)
    observed = explicit_packages(exported)
    missing = sorted(expected - observed)
    unexpected = sorted(observed - expected)
    if missing or unexpected:
        details = ['Workflow environment differs from the committed explicit lock.']
        if missing:
            details.append('Missing locked entries:\n  ' + '\n  '.join(missing))
        if unexpected:
            details.append('Unexpected exported entries:\n  ' + '\n  '.join(unexpected))
        details.append('The environment and exported evidence have been preserved; repair deliberately before retrying.')
        raise ValueError('\n'.join(details))
    return len(expected)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('lock', type=Path, help='Committed explicit package lock')
    parser.add_argument('exported', type=Path, help='Explicit export of the installed environment')
    args = parser.parse_args(argv)
    try:
        count = compare_environment(args.lock, args.exported)
    except (OSError, UnicodeError, ValueError) as error:
        print(f'Workflow environment check failed: {error}', file=sys.stderr)
        return 1
    print(f'Workflow environment matches the committed explicit lock ({count} packages).')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
