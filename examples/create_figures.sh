#!/bin/bash

here=$(pwd)
for examplefolder in `ls -d example-*`
do
    cd $examplefolder
    for cmguifile in `ls results/*.com`
    do
        # run cmgui script file and export screenshot
        echo $cmguifile
        cmgui $cmguifile
    done
    cd $here
done
