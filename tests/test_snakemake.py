"Tests snakemake execution via click CLI module."
import pytest
import tempfile
import shutil
import os
from pytest_dependency import depends

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
    conda_args = [] #["--conda-prefix", _tempdir]

    status = run_snakemake(conf, no_use_conda=True, verbose=True,
                           outdir=_tempdir, extra_args=[target] + conda_args + extra_args)
    return status

test_genomes = ['GCA_002691795.1_genomic.fna.gz', 'GCF_003143755.1_genomic.fna.gz']
test_proteomes = ['GB_GCA_002691795.1_protein.faa.gz', 'RS_GCF_003143755.1_protein.faa.gz']
nucl_databases = ['gtdb-nine.nucleotide-k21', 'gtdb-nine.nucleotide-k31', 'gtdb-nine.nucleotide-k51']
prot_databases = ['gtdb-nine.protein-k11', 'gtdb-nine.dayhoff-k19', 'gtdb-nine.hp-k33', 'gtdb-nine.hp-k42']

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


@pytest.mark.dependency(["test_nucleotide_sketch"])
@pytest.mark.parametrize("genome_file", test_genomes)
@pytest.mark.parametrize("database_name", nucl_databases)
def test_nucleotide_search_containment(request, genome_file, database_name):
    depends(request, [f"test_nucleotide_sketch[{g}]" for g in test_genomes])
    target = f"search/{genome_file}.x.{database_name}.matches.sig"
    status = _run_snakemake_test('config/test-nucl.yaml', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))

@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", test_genomes)
@pytest.mark.parametrize("database_name", prot_databases)
def test_translate_search_containment(request, genome_file, database_name):
    depends(request, [f"test_translate_sketch[{g}]" for g in test_genomes])
    target = f"search/{genome_file}.x.{database_name}.matches.sig"
    status = _run_snakemake_test('config/test-nucl.yaml', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", test_proteomes)
@pytest.mark.parametrize("database_name", prot_databases)
def test_protein_search_containment(request, genome_file, database_name):
    depends(request, [f"test_protein_sketch[{g}]" for g in test_proteomes])
    target = f"search/{genome_file}.x.{database_name}.matches.sig"
    status = _run_snakemake_test('config/test-prot.yaml', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


#@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", test_genomes)
@pytest.mark.parametrize("database_name", nucl_databases)
def test_nucleotide_contigs_taxonomy(request, genome_file, database_name):
    depends(request, [f"test_nucleotide_search_containment[{genome_file},{database_name}]"])
    target = f"classify/{genome_file}.x.{database_name}.contigs-tax.json"
    status = _run_snakemake_test('config/test-nucl.yaml', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", test_genomes)
@pytest.mark.parametrize("database_name", prot_databases)
def test_translate_contigs_taxonomy(request, genome_file, database_name):
    depends(request, [f"test_translate_search_containment[{genome_file}, {database_name}]"])
    target = f"classify/{genome_file}.x.{database_name}.contigs-tax.json"
    status = _run_snakemake_test('config/test-nucl.yaml', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", test_proteomes)
@pytest.mark.parametrize("database_name", prot_databases)
def test_protein_contigs_taxonomy(request, genome_file, database_name):
    depends(request, [f"test_protein_search_containment[{genome_file},{database_name}]"])
    target = f"classify/{genome_file}.x.{database_name}.contigs-tax.json"
    status = _run_snakemake_test('config/test-prot.yaml', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


