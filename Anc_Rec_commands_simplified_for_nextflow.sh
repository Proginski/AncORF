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

echo $query

# Append the neighbor short name to the toali file
#awk -v focal=$focal_species '$0 == focal {print ">"substr($0,1,5)FNR}' $species_list | head -c 10 > ${process_dir}/$out_seq
awk -F"," -v focal=$focal_species '$1 == focal {print ">"$2}' ${process_dir}/$phylip_names > ${process_dir}/$out_seq
awk -F"," -v focal=$focal_species '$1 == focal {print $2" CDS"}' ${process_dir}/$phylip_names > ${process_dir}/$names_n_types

# Then, first the focal (query) nucleotide sequence
grep $query -A1 GENOMES/${focal_species}_CDS.nfasta | grep -v \> >> ${process_dir}/${out_seq}

# Add the focal species in the names for the phylogenetic tree
#awk -v focal=$focal_species '$0 == focal {print ">"substr($0,1,5)FNR}' $species_list | head -c 10 > ${process_dir}/${query}_names.txt

# Then the appropriate neighbor sequences
cat $neighbors | while read neighbor; do

	echo $neighbor

	# In the BLASTp output of focal vs neighbor :
        # Only keep one row per query ($1) base on minimal evalue ($3)
        # Only keep this row if it satisfies the evalue ($3) threshold and print the subject sequence ($2)
	CDS_sbjct=$(grep $query BLASTP/${focal_species}_blastp_vs_${neighbor}.out \
	| awk 'BEGIN{FS=OFS="\t"} {if (a[$1] == "" || a[$1]>$3) {a[$1]=$3; data[$1]=$0}} END{for (i in a) print data[i]}' \
	| awk -F"\t" '$1 !~ /#/ && $3 < 0.001 {print $2}')

	# If there is a CDS match
	if [ ! -z $CDS_sbjct ]; then

		echo $CDS_sbjct

                # Append the neighbor to the names for the phylogenetic tree
                #awk -v neighbor=$neighbor '$0 == neighbor {print ">"substr($0,1,5)FNR}' $species_list | head -c 10 >> ${process_dir}/${query}_names.txt

                # Append the neighbor short name to the toali file
                #awk -v neighbor=$neighbor '$0 == neighbor {print ">"substr($0,1,5)FNR}' $species_list | head -c 10 >> ${process_dir}/$out_seq
                awk -F"," -v neighbor=$neighbor '$1 == neighbor {print ">"$2}' ${process_dir}/$phylip_names >> ${process_dir}/$out_seq
                awk -F"," -v neighbor=$neighbor '$1 == neighbor {print $2" CDS"}' ${process_dir}/$phylip_names >> ${process_dir}/$names_n_types
                
		# Elongate the range
                # First if not already done, get an FASTA with elongated CDS
		if [ ! -s ${process_dir}/${neighbor}_CDS_elongated.nfasta ]; then
			echo "${process_dir}/${neighbor}_CDS_elongated.nfasta does not exists or is empty. Generating it with ORFget -elongate 100 ..."
			python /SCRIPTS/ORFget_v2.py -fna GENOMES/${neighbor}.fna -gff GENOMES/${neighbor}.gff -o ${process_dir}/${neighbor}_CDS -type nfasta -features_include CDS -elongate 100
		fi

		# Only add the sequence not the correpsonding header
		grep $CDS_sbjct -A1 ${process_dir}/${neighbor}_CDS_elongated.nfasta | grep -v ">" >> ${process_dir}/$out_seq


	else
		echo "NO CDS match, looking for IGR match..."
	        # In the tBLASTn output of focal vs neighbor :
	        # Only keep one row per query ($1) base on minimal evalue ($3)
        	# Only keep this row if it satisfies the evalue ($3) threshold and print the subject sequence ($2)
        	IGR_line=$(grep $query TBLASTN/${focal_species}_TRG_tblastn_vs_${neighbor}.out \
        	| awk 'BEGIN{FS=OFS="\t"} {if (a[$1] == "" || a[$1]>$3) {a[$1]=$3; data[$1]=$0}} END{for (i in a) print data[i]}' \
        	| awk 'BEGIN{FS=OFS="\t"} $1 !~ /#/ && $3 < 0.001')

		IGR_sbjct=$(echo "$IGR_line" | awk -F"\t" '{print $2}')

	        # If there is an IGR match
	        if [ ! -z $IGR_sbjct ]; then

        	        echo $IGR_sbjct

                         # Appnd the neighbor to the names for the phylogenetic tree
                        #awk -v neighbor=$neighbor '$0 == neighbor {print ">"substr($0,1,5)FNR}' $species_list | head -c 10 >> ${process_dir}/${query}_names.txt

                        # Append the neighbor short name to the toali file
                        #awk -v neighbor=$neighbor '$0 == neighbor {print ">"substr($0,1,5)FNR}' $species_list | head -c 10 >> ${process_dir}/$out_seq
                        awk -F"," -v neighbor=$neighbor '$1 == neighbor {print ">"$2}' ${process_dir}/$phylip_names >> ${process_dir}/$out_seq
                        awk -F"," -v neighbor=$neighbor '$1 == neighbor {print $2" IGR"}' ${process_dir}/$phylip_names >> ${process_dir}/$names_n_types
			
			# Remove comment lines
			# In the second field of the tblastn output (matching IGR name), replace the pattern "nonalphanumericcharacter digigits - digits " by "tabulation digits tabulation digits"  so we can acces the start and the stop of the IGR
			# Elongate the hit as proposed by Papadopoulos and print the result as a one-line GFF.
			elongated_hit=$(awk 'BEGIN{FS=OFS="\t"} $1 !~ /#/' <(echo "$IGR_line") \
			| awk 'BEGIN{FS=OFS="\t"} {$2=gensub( /(.*)[^a-zA-Z0-9]([0-9]+)-([0-9]+)/ , "\\1\t\\2\t\\3" , 1 , $2)}1' \
			| awk 'BEGIN{FS=OFS="\t"}

			        function abs(x){return ((x < 0.0) ? -x : x)}
			        {
			        hit_seq_size = (abs($11-$10)+1)/3
			        positions_to_M = ($7+1)*3
			        positions_to_STOP = ($6 - $8)*3
			        reduce_from_STOP = (hit_seq_size - $12)*3
			        if ($10 < $11) {
					strand = "+"
			                first_pos = $3 + $10 - positions_to_M
			                last_pos  = $3 + $11 - reduce_from_STOP + positions_to_STOP
			        } else {
					strand = "-"
			                first_pos = $3 + $11 + reduce_from_STOP - positions_to_STOP
			                last_pos  = $3 + $10 + positions_to_M -1
			        }
				print $2,"ELONGATED_HIT","IGR",first_pos,last_pos,".",strand,".","X"
			}')


			# Get the corresponding sequence with ORFget and the option elongate

			# First if not already done, get a file with the chromosome lenghts
			if [ ! -s GENOMES/${neighbor}.genome ]; then faSize GENOMES/${neighbor}.fna -tab -detailed > GENOMES/${neighbor}.genome ; fi

			# Elongate the genomic range
			bedtools slop -i <(echo "$elongated_hit") -g GENOMES/${neighbor}.genome -b 100 > ${process_dir}/${query}_${neighbor}_elongated.gff
			# Get the corresponding sequence and append it to the alignment FASTA | wihtout the header (already written, see above)
			bedtools getfasta -fi GENOMES/${neighbor}.fna -bed ${process_dir}/${query}_${neighbor}_elongated.gff -s | grep -v ">" >> ${process_dir}/$out_seq


	        else
	                echo "NO IGR match."
	        fi
	fi
done

# Only retain cases where at least 3 sequences (1 focal + 2 neighbors) have been collected
seq_nb=$(cat ${process_dir}/$names_n_types | wc -l)

if [ $seq_nb -gt 2 ]; then
	echo "OK enough (${seq_nb}) sequences"
	names=$(awk '{print $1}' ${process_dir}/$names_n_types)
	echo "names = $names"
	# Generate the phylogenetic tree with only the species that have a match
	python /SCRIPTS/prune_tree.py -tree ${process_dir}/${tree} -names <(echo "$names") -out ${process_dir}/${query}.nwk
else
	echo "NOT enough (${seq_nb}) sequences"
	mv ${process_dir}/$out_seq ${process_dir}/fail
	mv ${process_dir}/$names_n_types ${process_dir}/fail_names_n_types
fi
