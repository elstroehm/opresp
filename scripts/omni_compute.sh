#!/bin/bash

#SBATCH --job-name=
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=16G
#SBATCH --mail-type=END
#SBATCH --mail-user=joel.stroehmann@student.uni-siegen.de
#SBATCH --array=0-26

LD_LIBRARY_PATH=/home/g039165/masterthesis/qsim/external/hdf5-1.14.1-2/hdf5/lib
export LD_LIBRARY_PATH
SECONDS=0
/home/g039165/masterthesis/qsim/build/g++/tls_response.out init/init_$SLURM_ARRAY_TASK_ID.json result/result_$SLURM_ARRAY_TASK_ID.h5
duration=$SECONDS
echo "$(($duration / 60)) Minuten und $(($duration % 60)) Sekunden Laufzeit."
