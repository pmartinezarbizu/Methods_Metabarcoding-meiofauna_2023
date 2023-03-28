#!/bin/bash
#
# copyright Pedro Martinez Arbizu & Sahar Khodami (2023)
# pmartinez@senckenberg.de
#  
### Trimm primers and adapters with bbmap
###use bbmap from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/

# run like: bash ./trimmadpat.sh  > log.txt 2>&1 &


#####################################
### EDIT THESE OPTIONS            ###
#####################################

	# bbmap path
	bbpath='/Dataexchange/bbmap'

	# raw sequence path
	RS="/Dataexchange/dada2batch/rawseq"
	
	
	#Primer V1V2
	primer='GCTTGTCTCAAAGATTAAGCC,GCCTGCTGCCTTCCTTGGA'

    # kmer size for searching primer
	# should be at most size of primer
	k=15
	
	# Hamming distance allowed
	hd=1

#####################################
### Do NOT EDIT BEXOND THIS LINE  ###
#####################################



	for f in $RS/*_R1_*.fastq.gz; do
	    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
 	    s=$(cut -d_ -f1 <<< "$f")
	   	   
	#for f in $RS/*_R1_*.fastq; do
	    #r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
 	    #s=$(cut -d_ -f1 <<< "$f")

	echo ====================================	
	echo "Step 1 : Trimm adapter overhang using paired-end reads and primer sequence"
	echo "processing $s"
	echo ====================================
	echo
	
	$bbpath/bbduk.sh in1=$f in2=$r out1=$s\_R1_tr.fastq out2=$s\_R2_tr.fastq literal=$primer \
	copyundefined ktrim=l k=$k hdist=$hd rcomp=t tbo overwrite=true
	
	echo 
	echo

	done

# move on with dada2
#	nohup Rscript --vanilla SGN_dada2_batch.r &

