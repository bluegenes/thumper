#ll:
#	@echo "You can run the following make targets: quicktest and test"

#quicktest:
	#py.test -k "not snakemake" tests

test:
	py.test tests/test_snakemake.py

verbosetest:
	py.test tests/test_snakemake.py -s
	#py.test tests
