default:
	make clean
	make run-tests-release
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

configure-example-debug=cmake -DOpenCMISSLibs_DIR=$(OPENCMISS_INSTALL_DIR) -DOPENCMISS_CONFIG_BUILD_TYPE=Debug ..

example-%:
	cd examples/$@ \
	&& mkdir -p debug \
	&& cd debug \
	&& cmake -DOpenCMISSLibs_DIR=$(OPENCMISS_INSTALL_DIR) -DOPENCMISS_CONFIG_BUILD_TYPE=Debug .. \
	&& make \
	&& make install \
	&& cd .. \
	&& bash run_example.sh debug \
	&& cd results \
	&& rm -f failed.tests results.summary \
	&& python compare_solutions.py

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
	rm -rf examples/example-0001/Makefile
	rm -rf examples/example-*/install_manifest.txt
	rm -rf examples/example-*/src/iron/example
	rm -rf examples/example-*/src/iron/CMakeFiles
	rm -rf examples/example-*/src/iron/Makefile
	rm -rf examples/example-*/src/iron/cmake_install.cmake
