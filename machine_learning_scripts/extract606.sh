#!/usr/bin/bash

###################################################################
################# Editing the HTML file into a fasta format #######
###################################################################
#First removing the 1st 14 lines, then,
#Removing html prefix of body from line 15 and retain the protein header

for filename in ~/Desktop/nehe/pat_data/*_fasta.mht
do
	echo "Working With File $filename"
	base=$(basename ${filename} _fasta.mht)
	echo "Base Name is $base"
	
	fasta=~/Desktop/nehe/pat_data/${base}.fa		#Destination of fasta file formats
	fasta2=~/Desktop/nehe/pat_data/raw/${base}_2.fa  
	sed '1,14d' $filename | sed 's/<BODY><PRE>//' | more | sed 's/<\/PRE><\/BODY><\/HTML>//' > $fasta 
	sed 's/&gt;/>/' $fasta > $fasta2 		#Replacing &gt with >
	mkdir ~/Desktop/nehe/pat_data/data_output/${base}  #Making directories where cellular protein files will be deposited.

echo "Done"
done


#######################################################################
############# Splitting the Multi-fasta files into single files #######
######################################################################
echo "Splitting the Multi-fasta files into single files"

for filename in ~/Desktop/nehe/pat_data/raw/*_2.fa
do
	echo "Working With File $filename"
        base=$(basename ${filename} _2.fa)
        echo "Base Name is $base"
	folder=~/Desktop/nehe/pat_data/data_output/${base}
	files=~/Desktop/nehe/pat_data/data_output/${base}/${base}.fa
	sed 's/ \r//' $filename > $files			#Removing spaces from the protein headers. Protein headers will be filenames
	echo "done"
	echo "FINISHED REMOVING SPACES FROM HEADERS!!!!"
done

echo "BEGINNING TO SPLIT THE MULTI-FASTA FILE"
wkdir=~/Desktop/nehe/pat_data/data_output/
cd $wkdir

for i in ./*/;
do 
	cd $i
	cat *.fa | awk '{
	if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
		print $0 > filename
	}'

echo "Done"
cd ../
done


