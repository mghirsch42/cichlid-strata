#!/bin/bash

n_bps=(0 2 4)
fit_type=("linear")
noise=(0.01 0.05)

n_samples=100
min_len=5

data_file="data/tropheini/100kb/LG5/Pseudosimochromis_marginatus.txt"

for dbp in "${n_bps[@]}"; do
    for dft in "${fit_type[@]}"; do
        for dn in "${noise[@]}"; do
            save_path="data/simulation/${dft}/nbps_${dbp}/noise_${dn}/"
            # save_path="temp/"
            echo $save_path

            if [ ! -d "$save_path" ]; then
                echo "Making save path ${save_path}"
                mkdir -p $save_path
            fi
            
            Rscript generate_simulated_data.R \
                --data_file $data_file \
                --save_path $save_path \
                --fit_type $dft \
                --n_bps $dbp \
                --min_len $min_len \
                --noise $dn \
                --n_samples $n_samples
        done
    done
done
