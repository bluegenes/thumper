from setuptools import setup, find_packages

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name = 'thumper-bio',
    version = "0.1",
    description="a tool for classifying genomes",
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
    package_data = { "thumper": ["thumper.snakefile", "*.yaml", "*.ipynb", ".yml"] },
    setup_requires = [ "setuptools>=38.6.0",
                       'setuptools_scm', 'setuptools_scm_git_archive' ],
    use_scm_version = {"write_to": "thumper/version.py"},
    install_requires = ['snakemake>=5.24', 'click>=7']
)
