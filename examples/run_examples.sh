#!/bin/bash -x

echo "This script will be used for testing / triggering Jenkins / showing a test summary from Jenkins test"

echo "Should also test serial / parallel execution!"

here=$(pwd)
target=$1
for examplefolder in `ls -d example-*`
do
    cd $examplefolder
    bash run_example.sh $target
    cd $here
done
pwd
