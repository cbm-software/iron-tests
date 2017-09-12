#!/bin/bash

for cmguifile in `ls results/*.com`
do
    # run cmgui script file and export screenshot
    echo $cmguifile
    cmgui $cmguifile
done
