.PHONY: build build-check cli-docs cli-docs-check fast test paml lint coverage benchmark check

build:
	Rscript tools/build_radte.R

build-check:
	Rscript tools/build_radte.R --check

cli-docs:
	Rscript tools/generate_cli_reference.R

cli-docs-check:
	Rscript tools/generate_cli_reference.R --check

fast: build-check cli-docs-check
	RADTE_TEST_PROFILE=fast Rscript test/run_tests.R

test: build-check cli-docs-check
	RADTE_TEST_PROFILE=full Rscript test/run_tests.R

paml: build-check
	RADTE_TEST_PROFILE=full RADTE_RUN_PAML_TESTS=true Rscript test/run_tests.R

lint:
	Rscript tools/lint.R

coverage:
	Rscript tools/coverage.R

benchmark:
	Rscript benchmark/benchmark_scaling.R

check: build-check cli-docs-check lint fast
