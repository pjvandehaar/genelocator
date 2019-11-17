#!/usr/bin/env python3
# to install: `pip install -e .`
# to install latest from pypi: `pip3 install -U --upgrade-strategy=eager --no-cache-dir genelocator`
# to upload to pypi: `./setup.py publish`
# to update dependencies: `pip3 install kpa && kpa pip-find-updates`, edit, `pip3 install -U --upgrade-strategy=eager genelocator`, test

from setuptools import setup, Command  # type: ignore
from pathlib import Path
import sys


def get_version():
    d = {}
    version_filepath = Path(__file__).absolute().with_name('genelocator') / '__version__.py'
    exec(version_filepath.read_text(), d)
    return d['version']


class PublishCommand(Command):
    # TODO: add option --init to skip querying pypi for a previous version
    description = 'Build and publish the package'
    user_options = [('initialize=', None, 'Ignore checking pypi version because package is new')]

    def initialize_options(self):
        self.initialize = None

    def finalize_options(self):
        pass

    @staticmethod
    def print_bold(s):
        print('\033[1m{0}\033[0m'.format(s))

    @classmethod
    def die(class_, s, exit_status=1):
        class_.print_bold(s)
        sys.exit(exit_status)

    @classmethod
    def get_return_code(class_, command, expected_codes=(0, 1)):
        import subprocess
        code = subprocess.run(command.split() if isinstance(command, str) else command).returncode
        if code not in expected_codes:
            class_.die("The command {} returned exit status {} but should have returned one of {}".format(command, code, expected_codes))
        return code

    def run(self):
        import json
        import subprocess
        from urllib.request import urlopen

        # die if there are unstaged changess
        if self.get_return_code('git diff-files --quiet') == 1:
            self.die('=> git workdir has changes => please either revert or stage them')
        # increment the version if it's the same as the PyPI version
        if not self.initialize:
            pypi_url = 'https://pypi.python.org/pypi/GeneLocator/json'
            latest_version = json.loads(urlopen(pypi_url).read())['info']['version']
            # Note: it takes pypi a minute to update the API, so this can be wrong.
            version = get_version()
            if latest_version == version:
                new_version_parts = version.split('.')
                new_version_parts[2] = str(1 + int(new_version_parts[2]))
                new_version = '.'.join(new_version_parts)
                self.print_bold('=> incrementing version {} -> {}'.format(version, new_version))
                Path('genelocator/__version__.py').write_text("version = '{}'\n".format(new_version))
                subprocess.run(['git', 'stage', 'genelocator/__version__.py'])
                version = new_version
        # commit any staged changes
        git_index_returncode = self.get_return_code('git diff-index --quiet --cached HEAD')
        if git_index_returncode == 1:
            self.print_bold('=> committing staged changes')
            subprocess.run(['git', 'commit', '-m', version])
        # delete ./dist/GeneLocator-* and repopulate it and upload to PyPI
        if Path('dist').exists() and list(Path('dist').iterdir()):
            setuppy = Path('setup.py').absolute()  # check that we are where we think we are before unlinking
            assert setuppy.is_file() and '47383945' in setuppy.read_text()
            for child in Path('dist').absolute().iterdir():
                if not child.name.startswith('GeneLocator-'):
                    self.die("Name of file inside dist/ doesn't begin with GeneLocator-, indicating that we might be in the wrong directory")
                self.print_bold('=> unlinking {}'.format(child))
                child.unlink()
        subprocess.run('python3 setup.py sdist bdist_wheel --universal'.split(), check=True)
        if not Path('~/.pypirc').expanduser().exists():
            self.print_bold('=> warning: you need a ~/.pypirc')
        try:
            subprocess.run('twine --version'.split())
        except FileNotFoundError:
            self.die('=> Run `pip3 install twine` and try again')
        subprocess.run('twine upload dist/*'.split(), check=True)
        if git_index_returncode == 1:
            self.print_bold('=> Now do `git push`.')
        sys.exit(0)


setup(
    name='GeneLocator',
    version=get_version(),
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
    cmdclass={
        'pub': PublishCommand
    },
)
