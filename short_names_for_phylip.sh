#!/bin/bash

while [ $# -gt 0 ]; do

   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        # echo $1 $2 // Optional to see the parameter:value result
   fi

  shift
done


# This script takes as input a file with a list of names/ID and returns a new file with a second column (separated by a comma)
# with the first 5 characters of each name/ID plus its rank in the file, il the limit of 9 total characters.
# It was intended to generate limited-size unique names for PHYLIP format.

mapfile -t arr < $in

for i in "${!arr[@]}"; do

	line=$(echo "${arr[i]:0:5}$((i+1))")
	echo "${arr[i]},${line:0:9}"

done > $out
