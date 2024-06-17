#!/bin/sh

Numbersimulation=100
dupcost_values=(2 5 10 20 50 70 100) 
losscost=1
threshold=70
dir_name="WGD1Simphy_S1_G100_Duprate-18"
#stats_out=$dir_name/stats_out_fix_loss_d100_l1_t70_GV2_1.csv
#gen_file_name="all_genetrees_edited.txt"
gen_file_name="applied_loss_fix_all_genetrees_edited.txt"
#gen_file_name="applied_loss_decBoth_all_genetrees_edited.txt"
#gen_file_name="applied_loss_decNum_all_genetrees_edited.txt"

str="D(lca-simphy),D(greedy-simphy),D(ultragreedy-simphy),Dif_square(simphy-LCA), Dif_square(simphy_greedy),C(simphy),DH(simphy),NBL(simphy),C(lca),DH(lca),NBL(lca),C(greedy),DH(greedy),NBL(greedy),C(ultragreedy),DH(ultragreedy),NBL(ultragreedy), max_nb_dup_species, sim number"
write_csv(){
    echo \"$1\",\"$2\",\"$3\",\"$4\",\"$5\",\"$6\",\"$7\",\"$8\",\"$9\",\"${10}\",\"${11}\",\"${12}\",\"${13}\",\"${14}\",\"${15}\",\"${16}\",\"${17}\",\"${18}\",\"${19}\" >> $stats_out
}

for dupcost in "${dupcost_values[@]}"
do
	stats_out=$dir_name/test_stats_out_d${dupcost}_l1_t70_GV2_1.csv
	echo $str > $stats_out
done

for (( i=1; i<=$Numbersimulation; i++ ))
do
	newname="sim_${i}"
	#maxheight=$(head -n 1 "$dir_name/$newname/maxheight.txt")
	maxheight="-"
	
	if [ -s $dir_name/$newname/s_tree.trees ]; then

		./segmentalreconcile_GV2_1 -d 2 -l 1 -gf $dir_name/$newname/$gen_file_name -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_simphy.txt -al simphy&
		./segmentalreconcile_GV2_1 -d 2 -l 1 -gf $dir_name/$newname/$gen_file_name -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_lca.txt -al lca&
		./segmentalreconcile_GV2_1 -d 2 -l 1 -gf $dir_name/$newname/$gen_file_name -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_greedy2.txt -al fastgreedy&
		./segmentalreconcile_GV2_1 -d 5 -l 1 -gf $dir_name/$newname/$gen_file_name -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_greedy5.txt -al fastgreedy&
		./segmentalreconcile_GV2_1 -d 10 -l 1 -gf $dir_name/$newname/$gen_file_name -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_greedy10.txt -al fastgreedy&
		./segmentalreconcile_GV2_1 -d 20 -l 1 -gf $dir_name/$newname/$gen_file_name -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_greedy20.txt -al fastgreedy&
		./segmentalreconcile_GV2_1 -d 50 -l 1 -gf $dir_name/$newname/$gen_file_name -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_greedy50.txt -al fastgreedy&
		./segmentalreconcile_GV2_1 -d 70 -l 1 -gf $dir_name/$newname/$gen_file_name -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_greedy70.txt -al fastgreedy&
		./segmentalreconcile_GV2_1 -d 100 -l 1 -gf $dir_name/$newname/$gen_file_name -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_greedy100.txt -al fastgreedy&
		
		wait
		
		
		
		file="$dir_name/$newname/out_simphy.txt" 
		j=1  
		while read -r line; do
			if [ "$j" -eq 2 ]; then
				simphy_cost2="$line"
			elif [ "$j" -eq 5 ]; then
				simphy_dupheight="$line"
			elif [ "$j" -eq 8 ]; then
				simphy_nblosses="$line"
			fi
			j=$((j+1))  
		done < $file  
		#echo "$simphy_cost, $simphy_dupheight, $simphy_nblosses"
		
		simphy_cost5=$((simphy_dupheight * 5)) 
		simphy_cost5=$((simphy_cost5 + simphy_nblosses)) 
		
		simphy_cost10=$((simphy_dupheight * 10)) 
		simphy_cost10=$((simphy_cost5 + simphy_nblosses)) 
		
		simphy_cost20=$((simphy_dupheight * 20)) 
		simphy_cost20=$((simphy_cost5 + simphy_nblosses)) 
		
		simphy_cost50=$((simphy_dupheight * 50)) 
		simphy_cost50=$((simphy_cost5 + simphy_nblosses)) 
		
		simphy_cost70=$((simphy_dupheight * 70)) 
		simphy_cost70=$((simphy_cost5 + simphy_nblosses)) 
		
		simphy_cost100=$((simphy_dupheight * 100)) 
		simphy_cost100=$((simphy_cost5 + simphy_nblosses)) 
		
		
		file="$dir_name/$newname/out_lca.txt" 
		j=1  
		while read -r line; do
			if [ "$j" -eq 2 ]; then
				lca_cost2="$line"
			elif [ "$j" -eq 5 ]; then
				lca_dupheight="$line"
			elif [ "$j" -eq 8 ]; then
				lca_nblosses="$line"
			fi
			j=$((j+1))  
		done < $file  
		#echo "$lca_cost, $lca_dupheight, $lca_nblosses"
		lca_cost5=$((lca_dupheight * 5)) 
		lca_cost5=$((lca_cost5 + lca_nblosses)) 
		
		lca_cost10=$((lca_dupheight * 10)) 
		lca_cost10=$((lca_cost5 + lca_nblosses)) 
		
		lca_cost20=$((lca_dupheight * 20)) 
		lca_cost20=$((lca_cost5 + lca_nblosses)) 
		
		lca_cost50=$((lca_dupheight * 50)) 
		lca_cost50=$((lca_cost5 + lca_nblosses)) 
		
		lca_cost70=$((lca_dupheight * 70)) 
		lca_cost70=$((lca_cost5 + lca_nblosses)) 
		
		lca_cost100=$((lca_dupheight * 100)) 
		lca_cost100=$((lca_cost5 + lca_nblosses)) 
		



		ultragreedy_nblosses="-"
		ultragreedy_dupheight="-"
		ultragreedy_cost="-"
		distance_ultragreedy_simphy="-"
		python compare_mapping.py $dir_name/$newname/out_lca.txt $dir_name/$newname/out_simphy.txt $dir_name/$newname/s_tree.newick $dir_name/$newname/comparison_lca_simphy.txt
		echo "comparison of LCA and simphy mapping is finished!"	
		distance_lca_simphy=$(head -n 22 "$dir_name/$newname/comparison_lca_simphy.txt" | tail -n 1)


		for dupcost in "${dupcost_values[@]}"
		do
			stats_out=$dir_name/test_stats_out_d${dupcost}_l1_t70_GV2_1.csv
			
			greedy_file="out_greedy${dupcost}.txt"
			comparison_greedy_simphy="comparison_greedy${dupcost}_simphy.txt"
			python compare_mapping.py $dir_name/$newname/$greedy_file $dir_name/$newname/out_simphy.txt $dir_name/$newname/s_tree.newick $dir_name/$newname/$comparison_greedy_simphy
			echo "comparison of greedy ${dupcost} and simphy mapping is finished!"	
			distance_greedy_simphy=$(head -n 22 "$dir_name/$newname/$comparison_greedy_simphy" | tail -n 1)


			file="$dir_name/$newname/out_greedy${dupcost}.txt" 
			j=1  
			greedy_cost=-1
			while read -r line; do
				if [ "$j" -eq 2 ]; then
					greedy_cost="$line"
				elif [ "$j" -eq 5 ]; then
					greedy_dupheight="$line"
				elif [ "$j" -eq 8 ]; then
					greedy_nblosses="$line"
				fi
				j=$((j+1))  
			done < $file  


			read SL SG Sl_Sq SG_Sq < <(bash count_nb_dup_to_calculate_diff_v2.sh $dir_name $newname $threshold $greedy_file)		

			if [ "$greedy_cost" -gt -1 ]; then
				simphy_cost=$(eval echo \$simphy_cost${dupcost})
				lca_cost=$(eval echo \$lca_cost${dupcost})
				write_csv $distance_lca_simphy $distance_greedy_simphy $distance_ultragreedy_simphy $Sl_Sq $SG_Sq $simphy_cost $simphy_dupheight $simphy_nblosses $lca_cost $lca_dupheight $lca_nblosses $greedy_cost $greedy_dupheight $greedy_nblosses $ultragreedy_cost $ultragreedy_dupheight $ultragreedy_nblosses $maxheight $i
			fi
			


		done
	else
		stats=-1
	fi
done


