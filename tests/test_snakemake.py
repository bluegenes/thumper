"Tests snakemake execution via click CLI module."
import pytest
import tempfile
import shutil
import os

from thumper.__main__ import run_snakemake
from . import pytest_utils as utils

# started from ctb testing implementation (spacegraphcats/test_snakemake.py)
#
# NOTE re dependencies (@pytest.mark.dependency):
# - These basically duplicate the snakemake dependencies.
# - they're there for convenience, because...
# - ...if these are wrong, the tests will still succeed, they just may
#   do some extra work in some tests & take longer.


def setup_module(m):
    global _tempdir
    _tempdir = tempfile.mkdtemp(prefix='thumper_test')

def teardown_module(m):
    global _tempdir
    try:
        shutil.rmtree(_tempdir, ignore_errors=True)
    except OSError:
        pass


def _run_snakemake_test(conf, target, extra_args=[]):
    conf = utils.test_file(conf)
    target = os.path.join(_tempdir, target)
    conda_args = ["--conda-prefix", _tempdir]

    status = run_snakemake(conf, no_use_conda=False, verbose=True,
                           outdir=_tempdir, extra_args=[target] + conda_args + extra_args)
    return status

test_genomes = ['2.fa.gz','63.fa.gz']
test_proteomes = ['GB_GCA_002691795.1_protein.faa.gz', 'RS_GCF_003143755.1_protein.faa.gz']

@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", test_genomes)
def test_nucleotide_sketch(genome_file):
    target = f"signatures/{genome_file}.nucleotide.sig"
    status = _run_snakemake_test("config/test-nucl.yaml", target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))

@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", test_genomes)
def test_translate_sketch(genome_file):
    target = f"signatures/{genome_file}.protein.sig"
    status = _run_snakemake_test("config/test-nucl.yaml", target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))

@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", test_proteomes)
def test_protein_sketch(genome_file):
    target = f"signatures/{genome_file}.protein.sig"
    status = _run_snakemake_test("config/test-prot.yaml", target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


#@pytest.mark.dependency(depends=['test_nucleotide_sketch'])
#def test_nucleotide_search_containment():
#    global _tempdir
#
#    prot_conf = utils.relative_file('tests/config/prot-test.yaml')
#    targets = ["output.test-prot/signatures/2.fa.gz.protein.sig",
#               "output.test-prot/signatures/63.fa.gz.protein.sig"]
#
#    status = test_run_snakemake(dory_conf, verbose=True, outdir=_tempdir,
#                           extra_args=targets)
#    assert status == 0
#    assert os.path.exists(os.path.join(_tempdir, target))


#@pytest.mark.dependency(depends=['test_translate_sketch'])
#def test_translate_search_containment():
#    global _tempdir
#
#    dory_conf = utils.relative_file('spacegraphcats/conf/dory-test.yaml')
#    target = 'dory_k21_r1/catlas.csv'
#    status = test_run_snakemake(dory_conf, verbose=True, outdir=_tempdir,
#                           extra_args=[target])
#    assert status == 0
#    assert os.path.exists(os.path.join(_tempdir, target))


#@pytest.mark.dependency(depends=['test_protein_sketch'])
#def test_protein_search_containment():
#    global _tempdir
#
#    dory_conf = utils.relative_file('spacegraphcats/conf/dory-test.yaml')
#    target = 'dory_k21_r1/catlas.csv'
#    status = test_run_snakemake(dory_conf, verbose=True, outdir=_tempdir,
#                           extra_args=[target])
#    assert status == 0
#    assert os.path.exists(os.path.join(_tempdir, target))

