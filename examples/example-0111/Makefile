default:
	make clean
	make release-fast
	make compare-fast
	make visualize-fast

fast:
	make

big:
	make clean
	make release-big
	make compare-big
	make visualize-big

debug-fast:
	bash run_example.sh debug fast

release-fast:
	bash run_example.sh release fast

debug-big:
	bash run_example.sh debug big

release-big:
	bash run_example.sh release big

compare:
	make compare-fast

compare-fast:
	cd results; rm -f failed.tests results.summary; matlab -nodesktop -nosplash -r compare_solutions

compare-big:
	echo ">>>ERROR: Automatic comparison not implemented"

visualize-fast:
	echo ">>>ERROR: Automatic visualization not implemented"

visualize-big:
	echo ">>>ERROR: Automatic visualization not implemented"

clean:
	rm -rf *.diag *.log *.aux *.out debug/ release/ results/current_run/*.exnode results/current_run/*.exelem doc/figures/*-eps-converted-to.pdf
