#!/bin/bash

# For a given CDS, this script looks for homologs among the CDS (BLASTp) and intergenic regions (i.e. IGR) (tBLASTn) of the neighbor species.
# It generates :
# - a FASTA file with the CDS (nucl) seq and the homologs (CDS or IGR)
# - a txt file with the names (species) and types (CDS or IGR) present in the FASTA file
# - a nwk file with evolutionnary relationships bewteen the selected species.


while [ $# -gt 0 ]; do

   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        # echo $1 $2 // Optional to see the parameter:value result
   fi

  shift
done

trap "exit 1" TERM
export TOP_PID=$$

cd $main_dir

echo $query


add_to_fasta () {

        if [ -z $header ]; then
                echo "*** NO HEADER for the sequence : $seq ***"
                kill -s TERM $TOP_PID
        else 
                if [ -z $seq ]; then
                        echo "*** EMPTY SEQUENCE for the header $header . This EXACT header might be absent of the file you are parsing ***"
                        kill -s TERM $TOP_PID
                else
                        echo $header >> ${process_dir}/${out_seq}
                        echo $seq >> ${process_dir}/${out_seq}
                fi
        fi

        header=""
        seq=""
}



# Append the focal (species of interest) short name to the output FASTA file and to the output txt file
header=$(awk -F"," -v focal=$focal_species '$1 == focal {print ">"$2}' ${process_dir}/$phylip_names)
awk -F"," -v focal=$focal_species '$1 == focal {print $2" CDS"}' ${process_dir}/$phylip_names > ${process_dir}/$names_n_types

# Add the focal (query) nucleotide sequence to the FASTA file
# First line : make sur to deal with a linearized (long format) input FASTA
# Second line : get the sequence without the header
seq=$(cat GENOMES/${focal_species}_CDS.nfasta | awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}' \
| grep -x ">${query}" -A1  | grep -v \>)

# Safely append header and seq to the output FASTA
add_to_fasta

# Then search for homologs in the neighbor sequences for each neighbor species.
cat $species_list | grep -v -x $focal_species | while read neighbor; do

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

		# Append the neighbor short name to the output FASTA file and to the output txt file
                header=$(awk -F"," -v neighbor=$neighbor '$1 == neighbor {print ">"$2}' ${process_dir}/$phylip_names)
                awk -F"," -v neighbor=$neighbor '$1 == neighbor {print $2" CDS"}' ${process_dir}/$phylip_names >> ${process_dir}/$names_n_types
                
		# Elongate the range
                # First if not already done, get a FASTA with elongated CDS for the neighbor (better futur alignment).
                elongated_fna=${process_dir}/${neighbor}_${CDS_sbjct}_elongated.nfasta
                
		if [ ! -s $elongated_fna ]; then
			echo "$elongated_fna does not exists or is empty. Generating it with ORFget -elongate 100 ..."
			
			# Define a subset of the neighbor GFF file that only contain (the) chromosome(s) mentionning ${CDS_sbjct}
			grep ${CDS_sbjct} GENOMES/${neighbor}.gff | awk -F"\t" '{print $1}' | sort | uniq > ${process_dir}/${neighbor}_${CDS_sbjct}_chrs.txt
			awk 'BEGIN{FS=OFS="\t"} {if(FNR == NR){data[$0] = 1} else{ if($1 in data) {print $0}}}' ${process_dir}/${neighbor}_${CDS_sbjct}_chrs.txt GENOMES/${neighbor}.gff > ${process_dir}/${neighbor}_${CDS_sbjct}.gff
			
			python /SCRIPTS/ORFget_v2.py -fna GENOMES/${neighbor}.fna -gff ${process_dir}/${neighbor}_${CDS_sbjct}.gff -o ${process_dir}/${neighbor}_${CDS_sbjct} -type nfasta -features_include CDS -elongate 100
			
			# make sur to deal with a linearized (long format) input FASTA
			cat $elongated_fna | awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}' > ${elongated_fna}_tmp
			mv ${elongated_fna}_tmp $elongated_fna
		fi

		# Only add the sequence not the correpsonding header.
		# As /SCRIPTS/ORFget_v2.py add a ".._mRNA" suffix that may not be present in the blastp subject name, look also for ">${CDS_sbjct}_.*mRNA"
		seq=$(grep -x -E ">${CDS_sbjct}|>${CDS_sbjct}_.*mRNA" $elongated_fna -A1  | grep -v ">")

		# Safely append header and seq to the output FASTA
		add_to_fasta
		

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

			# Append the neighbor short name to the output FASTA file and to the output txt file
                        header=$(awk -F"," -v neighbor=$neighbor '$1 == neighbor {print ">"$2}' ${process_dir}/$phylip_names)
                        awk -F"," -v neighbor=$neighbor '$1 == neighbor {print $2" IGR"}' ${process_dir}/$phylip_names >> ${process_dir}/$names_n_types
			
			# Remove comment lines
			# In the second field of the tblastn output (matching IGR name), replace the pattern "nonalphanumericcharacter digits - digits " by "tabulation digits tabulation digits"  so we can acces the start and the stop of the IGR
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
			if [ ! -s GENOMES/${neighbor}.genome ]; then /packages/faSize GENOMES/${neighbor}.fna -tab -detailed > GENOMES/${neighbor}.genome ; fi

			# Elongate the genomic range
			bedtools slop -i <(echo "$elongated_hit") -g GENOMES/${neighbor}.genome -b 100 > ${process_dir}/${query}_${neighbor}_elongated.gff
			# Get the corresponding sequence and append it to the alignment FASTA | wihtout the header (already written, see above)
			# (bedtools getfasta does generates long lines FASTA files)
			seq=$(bedtools getfasta -fi GENOMES/${neighbor}.fna -bed ${process_dir}/${query}_${neighbor}_elongated.gff -s | grep -v ">")
			
			# Safely append header and seq to the output FASTA
			add_to_fasta


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
