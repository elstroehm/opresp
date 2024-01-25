#!/bin/bash

g_exec_path=~/masterthesis/rework_qsim/build/clang++-18/tls_response.out

cd $1/init
for fname in *; do
    rname=${fname/init/result};
    rname=${rname/json/h5};
    $g_exec_path $fname ../results/$rname
done


