#!/bin/bash  
#
# copyright Pedro Martinez Arbizu & Sahar Khodami (2023)
# pmartinez@senckenberg.de
#  

#this master script will call following additional scripts in sequences

#trim adapters
bash trimmadapt.V1V2.sh

# analysis with DADA2
Rscript --vanilla SGN_dada2_batch_V1V2.r 

# copy files to input order
cp *.txt /Dataexchange/blastn/input/

# blast and annotate sequences
bash SGNpipeline_dada2_V1V2.sh

