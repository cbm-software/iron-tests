default:
	make clean
	mkdir -p results/current_run
	make release-fast
	make visualize-fast

fast:
	make

big:
	make clean 
	make release-big
	make visualize-big

release: release-fast

debug: debug-fast

create-figures: visualize-fast

compare:
	echo "example has no compare target, because it only contains a convergence study and no real tests"

debug-fast:
	bash run_example.sh debug fast

release-fast:
	bash run_example.sh release fast

debug-big:
	bash run_example.sh debug big

release-big:
	bash run_example.sh release big

visualize-fast:
	cd src/scripts; python time_convergence.py 64 3.0; python time_convergence.py 64 1.0; cd ../../results/current_run; cp *.eps ../../doc/figures

visualize-big:
	cd src/scripts; python time_convergence.py 64 10.0; python time_convergence.py 64 3.0; python time_convergence.py 64 1.0; cd ../../results/current_run; cp *.eps ../../doc/figures

clean:
	rm -rf *.diag *.log *.aux *.out debug/ release/
	rm -rf results/current_run
	rm -rf results/failed.tests results/results.summary
	rm -f  failed.tests results.summary
	rm -rf doc/figures/*-eps-converted-to.pdf
	rm -rf Time*
	rm -rf tmp*
	rm -rf Monodomain*
