#!/bin/bash

while [ $# -gt 0 ]; do

   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        # echo $1 $2 // Optional to see the parameter:value result
   fi

  shift
done

cd $main_dir

out_tree="${process_dir}/${tree}_checked"
cp ${process_dir}/${tree} $out_tree

# Only get node names
# Get those found several times
cat $out_tree | sed "s/[(),;]/\n/g" | sed "s/:.*//g" | awk '$0 != ""' | sort | uniq -c \
| awk '$1 > 1' | while read line; do

	occ=$(echo $line | awk '{print $1}')
	name=$(echo $line | awk '{print $2}')

	echo "$name is found $occ times"

	if grep -q $name $species_list ; then
		echo "$name is in the input list of species... The tree must be corrected with unique taxon labels. Cannot pursue."
		exit 1

	# Add the rank of each occurrence at the each of the given occurence of the label
	else
		# If the node ends with a quote, keep it ONLY WORKS FOR SIMPLE QUOTES
		if [ "${name: -1}" == "'" ]; then
			for i in $(seq 1 $occ); do
                        	sed -i "s/${name:0:-1}/${name:0:-1}_${i}/${i}" $out_tree
			done
		else
			for i in $(seq 1 $occ); do
				sed -i "s/${name}/${name}_${i}/${i}" $out_tree
			done
		fi
	fi
done

