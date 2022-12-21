rm -f "max_covs.txt"

ls ANC_ORF/ANCESTORS/*.frags | while read frags
do

	awk 'BEGIN{FS=OFS="\t"}
	{
	focal[NR]=$1
	start[NR]=$7
	end[NR]=$8
	}
	END{
	for (ind1 in start){
		max_cov=-1
		for (ind2 in start ){
			if(ind1 != ind2){
				printf "rowA : "start[ind1]"-"end[ind1]"\nrowB : "start[ind2]"-"end[ind2]" "
				if(start[ind1] > start[ind2]){max_start=start[ind1]}
				else{max_start=start[ind2]}

				if(end[ind1] < end[ind2]){min_end=end[ind1]}
	                        else{min_end=end[ind2]}

	#			printf "max_start : "max_start", min_end : "min_end

				cov= min_end - max_start +1
				if(cov > 0){
					print "cov : "cov
					if(cov > max_cov){ max_cov = cov }
				}
				else{ print "no cov" }
			}
		}
		print "max cov : "max_cov
		print focal[ind1],max_cov >> "max_covs.txt"
	}
	}' $frags

done
