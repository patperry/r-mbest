
check:
	Rscript -e 'library("devtools"); if (require('doMC')) registerDoMC(4); test(".")'

.PHONY: check
