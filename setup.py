#!/usr/bin/env python3
# to install: `pip install -e .`
# to install latest from pypi: `pip3 install -U --upgrade-strategy=eager --no-cache-dir genelocator`
# to upload to pypi: `./setup.py publish`
# to update dependencies: `pip3 install kpa && kpa pip-find-updates`, edit, `pip3 install -U --upgrade-strategy=eager genelocator`, test

from setuptools import setup
from pathlib import Path
import sys

d = {}
version_filepath = Path(__file__).absolute().with_name('genelocator') / '__version__.py'
exec(version_filepath.read_text(), d)
version = d['version']

if sys.argv[-1] in ['publish', 'pub']:
    import json
    import subprocess
    from urllib.request import urlopen

    # make sure there's no unstaged changess
    git_workdir_returncode = subprocess.run('git diff-files --quiet'.split()).returncode
    assert git_workdir_returncode in [0, 1]
    if git_workdir_returncode == 1:
        print('=> git workdir has changes => please either revert or stage them')
        sys.exit(1)
    # if the local version is the same as the PyPI version, increment it
    pypi_url = 'https://pypi.python.org/pypi/GeneLocator/json'
    latest_version = json.loads(urlopen(pypi_url).read())['info']['version']
    # Note: it takes pypi a minute to update the API, so this can be wrong.
    if latest_version == version:
        new_version_parts = version.split('.')
        new_version_parts[2] = str(1 + int(new_version_parts[2]))
        new_version = '.'.join(new_version_parts)
        print('=> incrementing version {} -> {}'.format(version, new_version))
        Path('genelocator/__version__.py').write_text("version = '{}'\n".format(new_version))
        version = new_version
        subprocess.run(['git', 'stage', 'genelocator/__version__.py'])
    # commit any staged changes
    git_index_returncode = subprocess.run('git diff-index --quiet --cached HEAD'.split()).returncode
    assert git_index_returncode in [0, 1]
    if git_index_returncode == 1:
        print('=> git index has changes')
        subprocess.run(['git', 'commit', '-m', version])
    # make sure there's a ~/.pypirc
    if not Path('~/.pypirc').expanduser().exists():
        print('=> warning: you need a ~/.pypirc')
    # delete ./dist/GeneLocator-* and repopulate it and upload to PyPI
    if Path('dist').exists() and list(Path('dist').iterdir()):
        setuppy = Path('setup.py').absolute()  # check that we are where we think we are before unlinking
        assert setuppy.is_file() and '47383945' in setuppy.read_text()
        for child in Path('dist').absolute().iterdir():
            assert child.name.startswith('GeneLocator-'), child
            print('=> unlinking', child)
            child.unlink()
    subprocess.run('python3 setup.py sdist bdist_wheel --universal'.split(), check=True)
    try:
        subprocess.run('twine --version'.split())
    except FileNotFoundError:
        print('=> Run `pip3 install twine` and try again')
        sys.exit(1)
    subprocess.run('twine upload dist/*'.split(), check=True)
    if git_index_returncode == 1:
        print('=> Now do `git push`.')
    sys.exit(0)

setup(
    name='GeneLocator',
    version=version,
    description="A library for finding the nearest gene to a genomic location",
    long_description=Path(__file__).with_name('README.md').read_text(),
    long_description_content_type='text/markdown',  # Optional (see note above)
    author="Peter VandeHaar",
    author_email="pjvandehaar@gmail.com",
    url="https://github.com/pjvandehaar/genelocator",
    classifiers=[
        'Programming Language :: Python :: 3 :: Only',
        'Operating System :: Unix',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    license='MIT',
    packages=['genelocator'],
    entry_points={
        'console_scripts': [
            'gene-locator=bin.command_line:main',
            # TODO: create (and add script for) gene-downloader, which will prepare the environment the first time
        ]},
    include_package_data=True,
    zip_safe=False,
    python_requires=">=3.5",
    setup_requires=[
        'pytest-runner~=5.0',
    ],
    install_requires=[
        'intervaltree~=3.0',
    ],
    tests_require=[
        'pytest~=5.0',
    ],
)
