#! /bin/bash
mkdir -p $2
for FILE in $(ls -1 "$1")
do
  echo "Converting $FILE"
  # Extract filename e.g my/dir/myfile.ext -> myfile
  filename=$(basename $FILE)
  filename=${filename%.*}

  perl software/gfftools/genbank2gff.pl $1/$FILE > $2/${filename}.gff3
done

