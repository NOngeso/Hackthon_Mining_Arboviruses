#! /bin/bash
count=1
for FILE in $(ls -1 "$1")
do
  echo "${count} $(realpath $1/$FILE)">> $2
  (( count++ ))
done