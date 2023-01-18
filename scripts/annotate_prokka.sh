#!/bin/bash
for FILE in $(ls -1 "$1")
do
  echo "Annotating $FILE"
  # Extract filename e.g my/dir/myfile.ext -> myfile
  filename=$(basename $FILE)
  filename=${filename%.*}

  prokka --kingdom Orthornavirae --outdir $2/$filename --prefix $filename $1/$FILE
done