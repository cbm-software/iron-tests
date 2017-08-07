default:
	make clean
	mkdir -p results/current_run
	make release-fast
	make compare
	make visualize-fast

fast:
	make

big:
	make clean 
	make release-big
	make compare
	make visualize-big

release: release-fast

debug: debug-fast

create-figures: visualize-fast

debug-fast:
	bash run_example.sh debug fast

release-fast:
	bash run_example.sh release fast

debug-big:
	bash run_example.sh debug big

release-big:
	bash run_example.sh release big

compare:
	cd results; python compare_solutions.py

visualize-fast:
	cd results; bash create_animations.sh fast

visualize-big:
	cd results; bash create_animations.sh big

clean:
	rm -rf *.diag *.log *.aux *.out debug/ release/
	rm -rf results/current_run
	rm -rf results/failed.tests results/results.summary
	rm -f  failed.tests results.summary
	rm -rf doc/figures/*-eps-converted-to.pdf
	rm -rf Time*
	rm -rf tmp*
	rm -rf Monodomain*
