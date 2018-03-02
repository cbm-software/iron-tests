default:
	make clean
	make example-0001
	make example-0001-u
	make example-0002
	make example-0003
	make example-0004
	make example-0005
	make example-0011
	make example-0012
	make example-0013
	make example-0101
	make example-0102
	make example-0111
	make example-0112
	make example-0302-u
	make example-0401
	make example-0402
	make example-0404-c
	make collect-results
	make create-figures
	make latexmk-once

tex:
	xelatex tests
	bibtex  tests
	xelatex tests
	xelatex tests

latexmk-once:
	latexmk -pdf tests

latexmk:
	latexmk -pdf -pvc tests

run-tests-debug:
	cd examples && bash run_examples.sh debug

run-tests-release:
	cd examples && bash run_examples.sh release

collect-results:
	cd examples; bash collect_results.sh

create-figures:
	cd examples; bash create_figures.sh

example-%:
	cd examples/$@ \
	&& mkdir -p debug \
	&& rm -rf debug/* \
	&& cd debug \
	&& cmake -DOpenCMISSLibs_DIR=$(OPENCMISS_INSTALL_DIR) -DOPENCMISS_CONFIG_BUILD_TYPE=Debug -DIRON_MESH_READER_ROOT=$(IRON_MESH_READER_ROOT) .. \
	&& make \
	&& make install \
	&& cd .. \
	&& bash run_example.sh debug \
	&& cd results \
	&& rm -f failed.tests results.summary \
	&& python compare_solutions.py | grep -iE "Passed tests:|Testing 2D|Testing 3D"

clean:
	rm -rf *aux doc/*aux *.bbl *.blg *.fdb_latexmk  *.fls *.lof *.log *.lot *.out *.toc
	rm -rf examples/example-*/doc/figures/current_run_*eps
	rm -rf examples/example-*/doc/figures/*-eps-converted-to.pdf
	rm -rf examples/example-*/doc/*.aux
	rm -rf examples/example-*/*.out
	rm -rf examples/example-*/*.diag
	rm -rf examples/example-*/debug
	rm -rf examples/example-*/release
	rm -rf examples/example-*/results/current_run
	rm -rf examples/example-*/results/failed.tests
	rm -rf examples/example-*/results/results.summary
	rm -rf failed.tests results.summary
	rm -rf examples/example-*/src/cheart/gmon.out
	rm -rf examples/example-*/CMakeCache.txt
	rm -rf examples/example-*/cmake_install.cmake
	rm -rf examples/example-*/CMakeFiles
	rm -rf examples/example-*/Makefile
	rm -rf examples/example-*/bin/*
	rm -rf examples/example-*/install_manifest.txt
	rm -rf examples/example-*/src/iron/example
	rm -rf examples/example-*/src/iron/CMakeFiles
	rm -rf examples/example-*/src/iron/Makefile
	rm -rf examples/example-*/src/iron/cmake_install.cmake
