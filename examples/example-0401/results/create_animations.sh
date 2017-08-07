# Loop over subdirectories of current_run and call plot3d.py script with $kind parameter.
# If kind is fast, the script will only generate images of the first and last data files, otherwise it will create animations.
kind=$1


for directory in `ls current_run`
do
    echo $directory
    python2.7 ../src/scripts/plot3d.py current_run/$directory $kind
    mv *.eps ../doc/figures/ || true
    mv *.png ../doc/figures/ || true
    mv *.mp4 ../doc/figures/ || true  # do not fail if there are no animations created (like in 'fast' target)
done

