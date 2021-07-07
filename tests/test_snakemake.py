"Tests snakemake execution via click CLI module."
import pytest
import tempfile
import shutil
import os
from pytest_dependency import depends
import io
import sys

from thumper.__main__ import run_snakemake
from snakemake.io import expand
from . import pytest_utils as utils

# NOTE: @pytest.mark.dependency dependencies enable
#       tests to pick up from prior test endpoint rather
#       than running full snakemake workflow up until target


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
    conf = utils.get_testfile(conf)
    target = os.path.join(_tempdir, target)

    sys.stdout, old_out = io.StringIO(), sys.stdout
    sys.stderr, old_err = io.StringIO(), sys.stderr
    try:
        status = run_snakemake(conf, no_use_conda=True, verbose=True,
                               outdir=_tempdir,
                               extra_args=[target] + extra_args)
    finally:
        sys.stdout, new_out = old_out, sys.stdout
        sys.stderr, new_err = old_err, sys.stderr

    return status, new_out.getvalue(), new_err.getvalue()

samples = ['GCA_002691795', 'GCF_003143755']

@pytest.mark.dependency()
@pytest.mark.parametrize("sample", samples)
def test_nucl_sketch(sample):
    target = f"sigs/nucleotide/{sample}.nucleotide.sig"
    status, out, err = _run_snakemake_test("config/test-nucl.yaml", target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency()
@pytest.mark.parametrize("sample", samples)
def test_translate_sketch(sample):
    target = f"sigs/translate/{sample}.protein.sig"
    status, out, err = _run_snakemake_test("config/test-translate.yaml", target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency()
@pytest.mark.parametrize("sample", samples)
def test_prot_sketch(sample):
    target = f"sigs/protein/{sample}.protein.sig"
    status, out, err= _run_snakemake_test("config/test-prot.yaml", target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=expand("test_nucl_sketch[{sample}]", sample=samples))
def test_nucl_zip():
    target = f"sigs/test-nucl.nucleotide.nucleotide.queries.zip"
    status, out, err = _run_snakemake_test("config/test-nucl.yaml", target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=expand("test_translate_sketch[{sample}]", sample=samples))
def test_translate_zip():
    target = f"sigs/test-translate.translate.protein.queries.zip"
    status, out, err = _run_snakemake_test("config/test-translate.yaml", target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=expand("test_prot_sketch[{sample}]", sample=samples))
def test_prot_zip():
    target = f"sigs/test-prot.protein.protein.queries.zip"
    status, out, err = _run_snakemake_test("config/test-prot.yaml", target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=["test_nucl_zip"])
@pytest.mark.parametrize("sample", samples)
def test_nucl_gather(sample):
    target = f"gather/{sample}.nucleotide.nucleotide-k21.gather.csv"
    status, out, err = _run_snakemake_test("config/test-nucl.yaml", target)

    print(out)
    print(err)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))

@pytest.mark.dependency(depends=["test_translate_zip"])
@pytest.mark.parametrize("sample", samples)
def test_translate_gather(sample):
    target = f"gather/{sample}.translate.protein-k10.gather.csv"
    status, out, err = _run_snakemake_test("config/test-translate.yaml", target)

    print(out)
    print(err)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=["test_prot_zip"])
@pytest.mark.parametrize("sample", samples)
def test_prot_gather(sample):
    target = f"gather/{sample}.protein.protein-k10.gather.csv"
    status, out, err = _run_snakemake_test('config/test-prot.yaml', target)

    print(out)
    print(err)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=expand("test_nucl_gather[{sample}]", sample=samples))
def test_nucl_gather_pathlist():
    target = f"gather/test-nucl.nucleotide.nucleotide-k21.gather-pathlist.txt"
    status, out, err = _run_snakemake_test("config/test-nucl.yaml", target)

    print(out)
    print(err)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=expand("test_translate_gather[{sample}]", sample=samples))
def test_translate_gather_pathlist():
    target = f"gather/test-translate.translate.protein-k10.gather-pathlist.txt"
    status, out, err = _run_snakemake_test("config/test-translate.yaml", target)

    print(out)
    print(err)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=expand("test_prot_gather[{sample}]", sample=samples))
def test_prot_gather_pathlist():
    target = f"gather/test-prot.protein.protein-k10.gather-pathlist.txt"
    status, out, err = _run_snakemake_test("config/test-prot.yaml", target)

    print(out)
    print(err)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends = ["test_nucl_gather_pathlist"])
def test_nucl_classify():
    target = f"classify/test-nucl.nucleotide.nucleotide-k21.classifications.csv"
    status, out, err = _run_snakemake_test('config/test-nucl.yaml', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends = ["test_translate_gather_pathlist"])
def test_translate_classify():
    target = f"classify/test-translate.translate.protein-k10.classifications.csv"
    status, out, err = _run_snakemake_test('config/test-translate.yaml', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends = ["test_prot_gather_pathlist"])
def test_protein_classify():
    target = f"classify/test-prot.protein.protein-k10.classifications.csv"
    status, out, err = _run_snakemake_test('config/test-prot.yaml', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))

