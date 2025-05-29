#!/bin/bash

data_bps=(0 2 4)
data_fit_type=("linear" "exponential" "quadratic")
data_noise=(0.01 0.05)

est_fit_type=("linear" "exponential" "quadratic")
tests=("t" "f" "both" "either")
min_bps=0
max_bps=20
min_len=5

for dbp in "${data_bps[@]}"; do
    for dft in "${data_fit_type[@]}"; do
        for dn in "${data_noise[@]}"; do
            data_path="data/simulation10/${dft}/nbps_${dbp}/noise_${dn}/"
            for eft in "${est_fit_type[@]}"; do
                for test in "${tests[@]}"; do
                    save_file="results/simulation10/data_${dft}/est_${eft}/nbps_${dbp}_noise_${dn}_test_${test}.csv"
                    Rscript run_simulated_loop.R \
                        --data_path $data_path \
                        --save_file $save_file \
                        --fit_type $eft \
                        --min_bps $min_bps \
                        --max_bps $max_bps \
                        --min_len $min_len \
                        --test $test \
                        --early_stop \
                        --quiet
                done
            done
        done
    done
done