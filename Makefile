RSCRIPT= Rscript --vanilla

all: dist

check:
	$(RSCRIPT) -e 'library("devtools"); if (require('doMC')) registerDoMC(4); test(".")'

clean:
	$(RSCRIPT) -e 'devtools::clean_dll(".")'

dist:
	mkdir -p dist && cd dist && R CMD build ..

install:
	$(RSCRIPT) -e 'devtools::install(".")'

.PHONY: all clean check dist install

.PHONY: check
