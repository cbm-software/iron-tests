default:
	make clean
	mkdir -p results/current_run
	make release
	make compare
	make create-figures

debug:
	bash run_example.sh debug

release:
	bash run_example.sh release

compare:
	cd results; rm -f failed.tests results.summary; python compare_solutions.py

create-figures:
	bash create_figures.sh

clean:
	rm -rf *.diag *.log *.aux *.out debug/ release/
	rm -rf results/current_run
	rm -rf results/failed.tests results/results.summary
	rm -f  failed.tests results.summary
	rm -rf doc/figures/*-eps-converted-to.pdf
