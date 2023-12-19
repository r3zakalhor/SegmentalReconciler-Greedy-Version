#!/bin/sh

Numbersimulation=100

dir_name="s9_100"
mkdir $dir_name
dupcost=5
losscost=0.5

stats_out=$dir_name/stats_out2.csv

str="D(lca-simphy),D(greedy-simphy),D(ultragreedy-simphy),C(simphy),DH(simphy),NBL(simphy),C(lca),DH(lca),NBL(lca),C(greedy),DH(greedy),NBL(greedy),C(ultragreedy),DH(ultragreedy),NBL(ultragreedy), max_nb_dup_species, sim number"
echo $str > $stats_out

write_csv(){
    echo \"$1\",\"$2\",\"$3\",\"$4\",\"$5\",\"$6\",\"$7\",\"$8\",\"$9\",\"${10}\",\"${11}\",\"${12}\",\"${13}\",\"${14}\",\"${15}\",\"${16}\",\"${17}\" >> $stats_out
}

for (( i=46; i<=$Numbersimulation; i++ ))
do
	./simphy -sb f:0.000001 -gb u:-20,-15 -gt f:0.0 -gg f:0.0 -gp f:0.0 -ld ln:gb,3 -lb f:ld -lt f:0.0 -lg f:0.0 -rs 1 -rl f:99 -rg 1 -o $dir_name -sp f:2 -su f:0.00001 -sg f:1 -sl U:20,50 -st f:1000000 -om 1 -v 3 -od 1 -op 1 -oc 1 -on 1 -ol 1 -cs $RANDOM
	cd $dir_name
	newname="sim_${i}"
	rm -r $newname
	mv 1 $newname
	cp ${dir_name}.command $newname/${newname}.command
	cd $newname 
	cat g_trees* >> all_genetrees.txt
	echo "gene trees are combined!"
	#bash ../../counter_dup_per_species.sh
	
	
	maxheight=0
	
	if [ -s s_tree.trees ]; then
		
		cd ..
		cd ..

		python post-order-labeling.py $dir_name/$newname/s_tree.trees $dir_name/$newname/s_tree.newick

		python map_gene_trees.py $dir_name/$newname/all_genetrees.txt $dir_name/$newname/all_genetrees_edited.txt $dir_name/$newname

		./segmentalreconcile_new3 -d $dupcost -l $losscost -gf $dir_name/$newname/all_genetrees_edited.txt -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_simphy.txt -al simphy
		file="$dir_name/$newname/out_simphy.txt" 
		j=1  
		while read -r line; do
			if [ "$j" -eq 2 ]; then
				simphy_cost="$line"
			elif [ "$j" -eq 5 ]; then
				simphy_dupheight="$line"
			elif [ "$j" -eq 8 ]; then
				simphy_nblosses="$line"
			fi
			j=$((j+1))  
		done < $file  
		#echo "$simphy_cost, $simphy_dupheight, $simphy_nblosses"
		
		
		./segmentalreconcile_new3 -d $dupcost -l $losscost -gf $dir_name/$newname/all_genetrees_edited.txt -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_lca.txt -al lca
		file="$dir_name/$newname/out_lca.txt" 
		j=1  
		while read -r line; do
			if [ "$j" -eq 2 ]; then
				lca_cost="$line"
			elif [ "$j" -eq 5 ]; then
				lca_dupheight="$line"
			elif [ "$j" -eq 8 ]; then
				lca_nblosses="$line"
			fi
			j=$((j+1))  
		done < $file  
		#echo "$lca_cost, $lca_dupheight, $lca_nblosses"
		
		./segmentalreconcile_new3 -d $dupcost -l $losscost -gf $dir_name/$newname/all_genetrees_edited.txt -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_greedy.txt -al greedy
		file="$dir_name/$newname/out_greedy.txt" 
		j=1  
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
		#echo "$greedy_cost, $greedy_dupheight, $greedy_nblosses"
		
		./segmentalreconcile_new3 -d $dupcost -l $losscost -gf $dir_name/$newname/all_genetrees_edited.txt -sf $dir_name/$newname/s_tree.newick -o $dir_name/$newname/out_ultragreedy.txt -al ultragreedy
		file="$dir_name/$newname/out_ultragreedy.txt" 
		j=1  
		while read -r line; do
			if [ "$j" -eq 2 ]; then
				ultragreedy_cost="$line"
			elif [ "$j" -eq 5 ]; then
				ultragreedy_dupheight="$line"
			elif [ "$j" -eq 8 ]; then
				ultragreedy_nblosses="$line"
			fi
			j=$((j+1))  
		done < $file  
		#echo "$ultragreedy_cost, $ultragreedy_dupheight, $ultragreedy_nblosses"

		python compare_mapping.py $dir_name/$newname/out_lca.txt $dir_name/$newname/out_simphy.txt $dir_name/$newname/s_tree.newick $dir_name/$newname/comparison_lca_simphy.txt
		echo "comparison of LCA and simphy mapping is finished!"	
		distance_lca_simphy=$(head -n 22 "$dir_name/$newname/comparison_lca_simphy.txt" | tail -n 1)
		#echo "distance_lca_simphy: $distance_lca_simphy"
		
		python compare_mapping.py $dir_name/$newname/out_greedy.txt $dir_name/$newname/out_simphy.txt $dir_name/$newname/s_tree.newick $dir_name/$newname/comparison_greedy_simphy.txt
		echo "comparison of greedy and simphy mapping is finished!"	
		distance_greedy_simphy=$(head -n 22 "$dir_name/$newname/comparison_greedy_simphy.txt" | tail -n 1)
		#echo "distance_greedy_simphy: $distance_greedy_simphy"
		
		python compare_mapping.py $dir_name/$newname/out_ultragreedy.txt $dir_name/$newname/out_simphy.txt $dir_name/$newname/s_tree.newick $dir_name/$newname/comparison_ultragreedy_simphy.txt
		echo "comparison of ultragreedy and simphy mapping is finished!"	
		distance_ultragreedy_simphy=$(head -n 22 "$dir_name/$newname/comparison_ultragreedy_simphy.txt" | tail -n 1)
		#echo "distance_ultragreedy_simphy: $distance_ultragreedy_simphy"

		write_csv $distance_lca_simphy $distance_greedy_simphy $distance_ultragreedy_simphy $simphy_cost $simphy_dupheight $simphy_nblosses $lca_cost $lca_dupheight $lca_nblosses $greedy_cost $greedy_dupheight $greedy_nblosses $ultragreedy_cost $ultragreedy_dupheight $ultragreedy_nblosses $maxheight $i

	else

		stats=-1
		cd ..
		cd ..

	
	fi

done
