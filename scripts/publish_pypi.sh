#!/usr/bin/env bash
# MagCluster PyPI Release Helper Script
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." >/dev/null 2>&1 && pwd )"
cd "$DIR"

echo "=== Cleaning previous build artifacts ==="
rm -rf dist build src/*.egg-info

echo "=== Building source distribution (sdist) and binary wheel ==="
python -m build

echo "=== Verifying packages with twine check ==="
twine check dist/*

echo "=== Artifacts ready for release ==="
ls -lh dist/

echo ""
echo "Package SHA256 Checksum for Bioconda:"
shasum -a 256 dist/*.tar.gz

echo ""
echo "============================================================"
echo "To publish to TestPyPI (recommended first):"
echo "  twine upload --repository testpypi dist/*"
echo ""
echo "To publish to Production PyPI:"
echo "  twine upload dist/*"
echo "  (Username: __token__, Password: your pypi-xxx token)"
echo "============================================================"
