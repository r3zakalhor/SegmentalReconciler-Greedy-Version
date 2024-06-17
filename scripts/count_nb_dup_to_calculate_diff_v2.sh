#!/bin/sh
declare -i bool
declare -i nbgenetrees
declare -a duparr
nbgenetrees=100
bool=1


dir_name="$1"
threshold="$3"
#dupcost="$4"
arr_file_path[1]="out_simphy.txt"
arr_file_path[2]="out_lca.txt"
arr_file_path[3]="$4"
arr_file_path[4]="out_ultragreedy.txt"

start_processing_text="<DUPS_PER_SPECIES>"
end_processing_text="</DUPS_PER_SPECIES>"

start_processing_text2="<SPECIESTREE>"
end_processing_text2="</SPECIESTREE>"


cd $dir_name
newname="$2"
cd $newname 

str="species id, Simphy dups, lca dups, greedy dups, dup in diff gene trees"




write_csv(){
    echo \"$1\",\"$2\",\"$3\",\"$4\",\"$5\"  >> $stats_out
}

    name="nb_dups.csv"
    stats_out=$name
    echo $str > $stats_out


	
	# Check if the file exists
	if [ -e "${arr_file_path[1]}" ]; then

		

		# Use a while loop to read the file line by line
		while IFS= read -r line; do

			if [[ "$line" == *"$start_processing_text2"* ]]; then
				# Set the flag to start processing lines
				process_lines2=true
				continue
			fi

			if [[ "$line" == *"$end_processing_text2"* ]]; then
				# Set the flag to start processing lines
				process_lines2=false
				continue
			fi

			if [ "$process_lines2" == true ]; then
				num_internal=$(echo "$line" | awk -F '[^0-9]+' '{for(i=1;i<=NF;i++)if($i!="")n=$i}END{print n}')
				# Print the result
				#echo "Root of the tree ID: $num_internal"
				# Extract and count the numbers within single quotes
				num_leaves=$(echo "$line" | grep -o "'[0-9]\+'")
				num_leaves=$(echo "$num_leaves" | wc -l)
				# Print the result
				#echo "Count of numbers within single quotes: $num_leaves"
				total=$(($num_leaves+$num_internal))
				#echo "Total: $total"
				for ((j = 0; j < total; j++)); do
    					duparr1[$j]=0
						duparr_per_gene_tree[$j]=0
				done
			fi

	



			# Check if the line contains the start processing text
			if [[ "$line" == *"$start_processing_text"* ]]; then
				# Set the flag to start processing lines
				process_lines=true
				continue
			fi

			# Check if the line contains the end processing text
			if [[ "$line" == *"$end_processing_text"* ]]; then
				# Set the flag to stop processing lines
				process_lines=false
				break
			fi

			# Process the lines between start and end
			if [ "$process_lines" == true ]; then
				sp_id=$(echo "$line" | awk '{print $1}' | tr -d '[]')
				sp_id_old=$(echo "$line" | awk '{print $1}' | tr -d '[]')
				if echo "$sp_id" | grep -q "'"; then
  					#echo "String contains single quotes: $sp_id"
					sp_id="${sp_id//\'/}"
					sp_id=$(($sp_id+$num_internal))
					#echo "new String contains single quotes: $sp_id"
				fi
				
				word_count=$(echo "$line" | wc -w)
				nb_dups=$(( (word_count - 1) / 2 ))
				duparr1[$sp_id]=$nb_dups
				uniqe_count=$(echo "$line" | awk -F'[()]' '{for (i=2; i<=NF; i+=2) print $i}' | sort -u | wc -l)
				#echo "$uniqe_count"
				duparr_per_gene_tree[$sp_id]=$uniqe_count

				#echo "$sp_id_old ${duparr_per_gene_tree[$sp_id]}"
			fi
		done < "${arr_file_path[1]}"
	else
		echo "File not found: ${arr_file_path[1]}"
	fi
	
    #echo "1 simulation counted!"    

########################################################################################################################

	# Check if the file exists
	if [ -e "${arr_file_path[2]}" ]; then

		

		# Use a while loop to read the file line by line
		while IFS= read -r line; do

			if [[ "$line" == *"$start_processing_text2"* ]]; then
				# Set the flag to start processing lines
				process_lines2=true
				continue
			fi

			if [[ "$line" == *"$end_processing_text2"* ]]; then
				# Set the flag to start processing lines
				process_lines2=false
				continue
			fi

			if [ "$process_lines2" == true ]; then
				num_internal=$(echo "$line" | awk -F '[^0-9]+' '{for(i=1;i<=NF;i++)if($i!="")n=$i}END{print n}')
				# Print the result
				#echo "Root of the tree ID: $num_internal"
				# Extract and count the numbers within single quotes
				num_leaves=$(echo "$line" | grep -o "'[0-9]\+'")
				num_leaves=$(echo "$num_leaves" | wc -l)
				# Print the result
				#echo "Count of numbers within single quotes: $num_leaves"
				total=$(($num_leaves+$num_internal))
				#echo "Total: $total"
				for ((j = 0; j < total; j++)); do
    					duparr2[$j]=0
				done
			fi

	



			# Check if the line contains the start processing text
			if [[ "$line" == *"$start_processing_text"* ]]; then
				# Set the flag to start processing lines
				process_lines=true
				continue
			fi

			# Check if the line contains the end processing text
			if [[ "$line" == *"$end_processing_text"* ]]; then
				# Set the flag to stop processing lines
				process_lines=false
				break
			fi

			# Process the lines between start and end
			if [ "$process_lines" == true ]; then
				sp_id=$(echo "$line" | awk '{print $1}' | tr -d '[]')

				if echo "$sp_id" | grep -q "'"; then
  					#echo "String contains single quotes: $sp_id"
					sp_id="${sp_id//\'/}"
					sp_id=$(($sp_id+$num_internal))
					#echo "new String contains single quotes: $sp_id"
				fi
				
				word_count=$(echo "$line" | wc -w)
				nb_dups=$(( (word_count - 1) / 2 ))
				duparr2[$sp_id]=$nb_dups
				#echo "$sp_id ${duparr${i}[$sp_id]}"
			fi
		done < "${arr_file_path[2]}"
	else
		echo "File not found: ${arr_file_path[2]}"
	fi
	

    #echo "2 simulation counted!"   
	
	
	###################################################################################################################################
	
		# Check if the file exists
	if [ -e "${arr_file_path[3]}" ]; then

		

		# Use a while loop to read the file line by line
		while IFS= read -r line; do

			if [[ "$line" == *"$start_processing_text2"* ]]; then
				# Set the flag to start processing lines
				process_lines2=true
				continue
			fi

			if [[ "$line" == *"$end_processing_text2"* ]]; then
				# Set the flag to start processing lines
				process_lines2=false
				continue
			fi

			if [ "$process_lines2" == true ]; then
				num_internal=$(echo "$line" | awk -F '[^0-9]+' '{for(i=1;i<=NF;i++)if($i!="")n=$i}END{print n}')
				# Print the result
				#echo "Root of the tree ID: $num_internal"
				# Extract and count the numbers within single quotes
				num_leaves=$(echo "$line" | grep -o "'[0-9]\+'")
				num_leaves=$(echo "$num_leaves" | wc -l)
				# Print the result
				#echo "Count of numbers within single quotes: $num_leaves"
				total=$(($num_leaves+$num_internal))
				#echo "Total: $total"
				for ((j = 0; j < total; j++)); do
    					duparr3[$j]=0
				done
			fi

	



			# Check if the line contains the start processing text
			if [[ "$line" == *"$start_processing_text"* ]]; then
				# Set the flag to start processing lines
				process_lines=true
				continue
			fi

			# Check if the line contains the end processing text
			if [[ "$line" == *"$end_processing_text"* ]]; then
				# Set the flag to stop processing lines
				process_lines=false
				break
			fi

			# Process the lines between start and end
			if [ "$process_lines" == true ]; then
				sp_id=$(echo "$line" | awk '{print $1}' | tr -d '[]')

				if echo "$sp_id" | grep -q "'"; then
  					#echo "String contains single quotes: $sp_id"
					sp_id="${sp_id//\'/}"
					sp_id=$(($sp_id+$num_internal))
					#echo "new String contains single quotes: $sp_id"
				fi
				
				word_count=$(echo "$line" | wc -w)
				nb_dups=$(( (word_count - 1) / 2 ))
				duparr3[$sp_id]=$nb_dups
				#echo "$sp_id ${duparr${i}[$sp_id]}"
			fi
		done < "${arr_file_path[3]}"
	else
		echo "File not found: ${arr_file_path[3]}"
	fi
	
	
    #echo "3 simulation counted!"   
	
	diff_SL=0
	diff_SG=0
	diff_SL_Sq=0
	diff_SG_Sq=0

    for k in ${!duparr1[@]}; do
        write_csv $k ${duparr1[$k]} ${duparr2[$k]} ${duparr3[$k]} ${duparr_per_gene_tree[$k]}
		if (( ${duparr_per_gene_tree[$k]} >= threshold )); then
			SL=${duparr1[$k]}-${duparr2[$k]}
			SG=${duparr1[$k]}-${duparr3[$k]}
			diff_SL=$((diff_SL+SL))
			diff_SG=$((diff_SG+SG))
			diff_SL_Sq=$((diff_SL_Sq+$((SL**2))))
			diff_SG_Sq=$((diff_SG_Sq+$((SG**2))))
		fi
    done

diff_SL_Sq=$(echo "sqrt($diff_SL_Sq)" | bc -l)
diff_SG_Sq=$(echo "sqrt($diff_SG_Sq)" | bc -l)

echo "$diff_SL $diff_SG $diff_SL_Sq $diff_SG_Sq"

cd ..
cd ..