#!/bin/bash

if [ -d $1 ]; then
    echo "Simulation mit diesem Namen existiert berteits"
else
    mkdir $1;
    mkdir $1/plots;
    mkdir $1/init;
    mkdir $1/results;
fi
