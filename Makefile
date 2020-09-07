#ll:
#	@echo "You can run the following make targets: quicktest and test"

quicktest:
	pytest tests/test_search_utils.py
	#py.test -k "not snakemake" tests

snaketest:
	py.test tests/test_snakemake.py

snaketest_v:
	py.test tests/test_snakemake.py -s

snaketest_nucl:
	pytest tests/test_snakemake.py -v -s -k "002691795 and nucleotide"

snaketest_translate:
	pytest tests/test_snakemake.py -v -s -k "002691795 and nucleotide or translate"

snaketest_prot:
	pytest tests/test_snakemake.py -v -s -k "002691795 and protein"

#alltests:	
    #py.test tests
