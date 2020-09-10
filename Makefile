#all:
#	@echo "You can run the following make targets: quicktest and test"


quicktest:
	pytest -k "not snakemake" tests

snaketest:
	pytest tests/test_snakemake.py

snaketest_v:
	pytest tests/test_snakemake.py -s

snaketest_nucl:
	pytest tests/test_snakemake.py -v -s -k "002691795 and nucleotide"

snaketest_translate:
	pytest tests/test_snakemake.py -v -s -k "002691795 and nucleotide or translate"

snaketest_prot:
	pytest tests/test_snakemake.py -v -s -k "002691795 and protein"

alltests:
	pytest tests
