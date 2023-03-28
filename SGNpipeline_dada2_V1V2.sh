#!/bin/bash                                                                       
#
# copyright Pedro Martinez Arbizu & Sahar Khodami (2023)
# pmartinez@senckenberg.de
#
# SGN pipeline for mapping dada2 results to blast hits
# create 5 subfolders
# /scripts (holding the scripts, execute from this directory)
# /input (save the dada2 files here)
# /output for hoding the output files (note that they will be overwritten if present)
# /output/blast inside the output folder for holding the output of blast and the final Taxon table
# /blastdb for holding your blast databases and own library
#
#
#	edit variable THREADS = number of threads or cores for parallel computing
     

# run in ubuntu with bash ./SGNpipeline_dada2_V1V2.sh

# to log the output and errors to a file do:
#  bash ./SGNpipeline_dada2_V1V2.sh  >> ../output/logfile.txt 2>&1 &


#####################################
### EDIT THESE OPTIONS            ###
#####################################
# running options
	THREADS=60
	
# DADA2 files
	DADA2OUT=nonchim.dada2.txt # file with non chimeric sequences
	DADA2Tab=full.nonchim.dada2.txt # Community table from dada2 without col and row names

# Amplicon option
	GENE=V1V2 # Fragment name
	#GENE=COI
# options for blast
	NBALIGN=10 # how many hits to keep from genebank blast
	#OWNLIB=Ref_Lib_18S_270519.fasta  # the name of your own library concatenated with genebank blast hits
	#OWNLIB=Ref_Lib_18S_220920.fasta
	OWNLIB=Ref_Lib_18S_060421.fasta


	BLASTP='/usr/bin'
	BLASTDB='/Dataexchange/blastn'
#####################################
### DO NOT EDIT BEYOND THIS LINE  ###
#####################################


	
# cleaning up
	echo
	echo ====================================
	echo Removing old files from output directory
	echo ====================================

	date
	
	rm -rfv $BLASTDB/output
	mkdir $BLASTDB/output
	mkdir $BLASTDB/output/blast
	mkdir $BLASTDB/output/taxa
	
	echo
	echo ====================================
	echo Processing gene $GENE from DADA2
	echo ====================================

# Enter subdirectory with files                                         

	echo
	echo Checking FASTQ format version 
	cd $BLASTDB/input

# edit file to create fasta format
	sed -i 's/"x"//g' $DADA2OUT #exclude initial "x"
	sed -i '/^\s*$/d' $DADA2OUT #exclude empty line
	sed -i 's/"//g' $DADA2OUT # exclude all "
	sed -i 's/^/>/g' $DADA2OUT # create a > at begin of line
	sed -E -i  's/[0-9]{1,10}/&\n/' $DADA2OUT # all new line after last numeral
	sed -i 's/^ //g' $DADA2OUT # remove space at begin of line

# edit community table file to remove colum names
	sed -E -i 's/[AGTC"]{1,1000}/&\n/' $DADA2Tab
	sed  -i '/^"/d' $DADA2Tab


# Blast against Genebank, combine with own library and blast again
# Change the variant name with Taxon name
 
	echo ###################################	
	echo query agains blastdb
	echo ###################################	
	
	
	cd $BLASTDB
	blastn -db nt -query $BLASTDB/input/$DADA2OUT \
	-out $BLASTDB/output/blast/all.otus.$GENE.dada2.blast.txt \
	-num_alignments $NBALIGN  -num_threads $THREADS \
	-outfmt "6  qseqid pident length  sblastnames sscinames sacc evalue sseq"


	
	echo ###################################	
	echo deduplicating fasta from blast results
	echo ###################################	

	awk -F'\t' '{print ">"$6"|"$4"|"$5}'  \
	$BLASTDB/output/blast/all.otus.$GENE.dada2.blast.txt \
	| sort -uk1,1  > $BLASTDB/output/blast/all.otus.$GENE.dada2.blast.dedup.txt

	# exclude hit with low taxonomic resolution
	grep -v  'eukaryotes|uncultured\|eukaryotes|eukaryote\|animals|uncultured\|fungi|uncultured\|NA|uncultured_organism\|eukaryotes|marine_zooplankton_environmental_sample\|eukaryotes|marine_copepoda_environmental_sample\|eukaryotes|copepoda_environmental_sample' \
	$BLASTDB/output/blast/all.otus.$GENE.dada2.blast.dedup.txt > \
	$BLASTDB/output/blast/all.otus.$GENE.dada2.blast.dedup.clean.txt



	# create list of accession numbers
	sed 's/ /_/g' $BLASTDB/output/blast/all.otus.$GENE.dada2.blast.dedup.clean.txt \
	| awk -F'[>\|]' '{print $2}' > $BLASTDB/output/blast/all.otus.$GENE.dada2.blast.sacc.txt	 

	
    # retrieve full sequences with accession numbers
	blastdbcmd -db nt \
	-entry_batch $BLASTDB/output/blast/all.otus.$GENE.dada2.blast.sacc.txt \
	-outfmt "%a|%s"  > $BLASTDB/output/blast/all.otus.$GENE.dada2.blast.fullseq.txt	

	# remove version number
	sed -i -e 's/\..\?//g' $BLASTDB/output/blast/all.otus.$GENE.dada2.blast.fullseq.txt


    echo ###################################	
	echo deduplicate and create fasta_file from blast results
	echo ###################################	
	
	awk -F'[>\|]' 'NR==FNR{a[$1]=$2;next} ($2 in a) {print ">"$2"|"$3"|"$4"\n"a[$2]}'  \
	$BLASTDB/output/blast/all.otus.$GENE.dada2.blast.fullseq.txt \
	$BLASTDB/output/blast/all.otus.$GENE.dada2.blast.dedup.clean.txt \
	> $BLASTDB/output/blast/all.otus.$GENE.dada2.blast.dedup.fasta	
	
	#change space to _	
	sed -i 's/ /_/g' $BLASTDB/output/blast/all.otus.$GENE.dada2.blast.dedup.fasta  	
	
	echo ###################################
	echo concatenate blast hits with own library
	echo ###################################

	cat $BLASTDB/output/blast/all.otus.$GENE.dada2.blast.dedup.fasta $OWNLIB \
	    > finaldb.fasta

	echo ###################################
	echo making blast database ...
	echo ###################################

	makeblastdb -in finaldb.fasta -dbtype nucl


	echo ###################################
	echo quey against own blast db all.otus.$GENE.dada2.fasta ...
	echo ###################################

	blastn -db finaldb.fasta -query \
	    $BLASTDB/input/$DADA2OUT \
	    -out $BLASTDB/output/blast/all.otus.$GENE.dada2.final.blast.txt -num_alignments 1  \
	    -num_threads $THREADS \
	    -outfmt "10  qseqid pident qcovs evalue length sacc"

	#replace "|" ,  "," and ";" by space
	sed -i 's/|/ /g;s/,/ /g;s/;/ /g' $BLASTDB/output/blast/all.otus.$GENE.dada2.final.blast.txt 


	echo ###################################	
	echo deduplicating blast results extracting only first species name provided by blast
	echo ###################################	
	
	awk '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8}' $BLASTDB/output/blast/all.otus.$GENE.dada2.final.blast.txt | \
	sort -uk1,1 -V > $BLASTDB/output/blast/all.otus.$GENE.dada2.final.dedup.blast.txt

	#format as table
	column -t $BLASTDB/output/blast/all.otus.$GENE.dada2.final.dedup.blast.txt \
	    > $BLASTDB/output/blast/all.otus.$GENE.dada2.final.blast.tab.txt

	column -t $BLASTDB/input/$DADA2Tab \
	    > $BLASTDB/output/blast/all.otutab.$GENE.dada2.tab.txt


	#create fasta file with original OTU sequence and description from blast
	grep '^[AGTC]' $BLASTDB/input/$DADA2OUT > $BLASTDB/output/blast/otuseqs.$GENE.dada2.txt

	#remove ASV which had "No hit"
	Rscript --vanilla $BLASTDB/scripts/removeR.r


	paste -d '\t' $BLASTDB/output/blast/all.otus.$GENE.dada2.final.dedup.blast.txt $BLASTDB/output/blast/otuseqnew.txt \
	| awk  '{print ">"$1"|"$6"|"$7"|"$8"|"$1"|"$2"|"$3"|"$4"|"$5"\n"$9}' \
	>  $BLASTDB/output/taxa/All.OTUS.blast.$GENE.dada2.fasta


	#extract by taxon
	declare -a taxa=($(awk -F '|' '{print $3}' $BLASTDB/output/taxa/All.OTUS.blast.$GENE.dada2.fasta | sort -uk1,1))
	
	for taxon in "${taxa[@]}"; do
	grep -A1 --no-group-separator $taxon $BLASTDB/output/taxa/All.OTUS.blast.$GENE.dada2.fasta \
	 > $BLASTDB/output/taxa/$taxon.$GENE.dada2.fasta

	done #done extract by taxon

	# clean up temporary files
	#rm ../output/blast/otuseqs.$GENE.dada2.txt
	#rm ../output/blast/all.otus.$GENE.dada2.final.blast.tab.txt 
	#rm ../output/blast/all.otus.$GENE.dada2.final.dedup.blast.txt
	#rm ../output/blast/all.otus.$GENE.dada2.blast.dedup.fasta
	#rm ../output/blast/all.otutab.$GENE.dada2.tab.txt 
	#rm ../output/blast/all.otus.$GENE.dada2.blast.sacc.txt 
	#rm ../output/blast/all.otus.$GENE.dada2.blast.fullseq.txt	


#done # done OTUSIM loop


# Game over
	echo
	echo Done
	date
