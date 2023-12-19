#!/bin/sh
declare -i bool
declare -i nbgenetrees
declare -a duparr
nbgenetrees=100
bool=1

Numbersimulation=100

dir_name="s9_100"
file_path="out_simphy.txt"
out_path="nb_dup_simphy_reconciliation_21"

start_processing_text="<DUPS_PER_SPECIES>"
end_processing_text="</DUPS_PER_SPECIES>"

start_processing_text2="<SPECIESTREE>"
end_processing_text2="</SPECIESTREE>"


cd $dir_name

mkdir $out_path  # Create the directory if it doesn't exist

str="species id,nb dup"




write_csv(){
    echo \"$1\",\"$2\" >> $stats_out
}

for (( i=1; i<=$Numbersimulation; i++ )); do

    
    
    name="nb_dup_sim_${i}.csv"
    stats_out=$out_path/$name
    echo $str > $stats_out

    newname="sim_${i}"
    cd $newname 

    #duparr=(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0) # size of 120


	
	# Check if the file exists
	if [ -e "$file_path" ]; then

		

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
    					duparr[$j]=0
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
				duparr[$sp_id]=$nb_dups
				#echo "$sp_id ${duparr[$sp_id]}"
			fi
		done < "$file_path"
	else
		echo "File not found: $file_path"
	fi
	
    cd ..
	
    echo "$i simulation counted!"    

    for k in ${!duparr[@]}; do
        write_csv $k ${duparr[$k]}
    done


done
