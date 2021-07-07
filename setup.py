from setuptools import setup, find_packages

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name = 'thumper-bio',
    description="workflow for classifying genomes using sourmash",
    url="https://github.com/bluegenes/thumper",
    author="N. Tessa Pierce and C. Titus Brown",
    author_email="ntpierce@gmail.com, titus@idyll.org",
    license="BSD 3-clause",
    packages = find_packages(),
    classifiers = CLASSIFIERS,
    entry_points = {'console_scripts': [
        'thumper  = thumper.__main__:main'
        ]
    },
    include_package_data=True,
    package_data = { "thumper": ["Snakefile","*.ipynb", ".yml"] },
    setup_requires = [ 'setuptools>=38.6.0',
                       'setuptools_scm',
                       'setuptools_scm_git_archive',
                       'pytest-runner',
                       ],
    use_scm_version = {"write_to": "thumper/version.py"},
    install_requires = ['snakemake>=6.5.2', 'click>=7,<8', 'pandas>1,<2'],
    long_description=long_description,
    long_description_content_type="text/markdown",
)
