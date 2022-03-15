#!/bin/bash
#submit using: bsub < /home/raanreye/Submit_FastqFragment.sh
#BSUB -J Well_job
#BSUB -o Well_job.%J.out
#BSUB -e Well_job.%J.error
#BSUB -M 200000 # Reqesting 10GB RAM
#BSUB -R "span[hosts=1] rusage [mem=200000]"
#BSUB -q normal #input gpu if using gpu

module load python/3.9.1
python /project/shafferslab/Rodrigo/maester_scripts/Fastq_sequence_fragment_and_clean.py