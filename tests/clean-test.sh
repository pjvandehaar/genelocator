#!/bin/bash
set -euo pipefail
_readlinkf() { perl -MCwd -le 'print Cwd::abs_path shift' "$1"; }
script_dir="$(cd "$(dirname "$(_readlinkf "${BASH_SOURCE[0]}")")" && echo "$PWD")"
repo_dir="$(dirname "$script_dir")"

cd "$repo_dir"
echo "=> Making wheel"
if [[ -d GeneLocator.egg-info ]]; then rm -r GeneLocator.egg-info; fi
if [[ -d build/lib/genelocator ]]; then rm -r build; fi
if [[ -d dist/GeneLocator* ]]; then rm -r dist; fi
python3 setup.py sdist bdist_wheel --universal

cd /tmp
rm -rf ./clean-test-venv
echo "=> Making venv"
python3 -m venv clean-test-venv
./clean-test-venv/bin/pip3 install "$repo_dir"/dist/GeneLocator-*.whl
./clean-test-venv/bin/pip3 install pytest
echo "=> Running pytest"
./clean-test-venv/bin/pytest "$repo_dir"/tests/

echo "=> Success"
