#!/bin/sh

Numbersimulation=100

dir_name="WGD2Simphy_S2_G100"
gen_file_name="all_genetrees_edited.txt"
#gen_file_name="applied_loss_fix_all_genetrees_edited.txt"
#gen_file_name="applied_loss_decBoth_all_genetrees_edited.txt"
#gen_file_name="applied_loss_decNum_all_genetrees_edited.txt"


for (( i=1; i<=$Numbersimulation; i++ ))
do
	./simphy -rs 1 -rl f:1 -rg 1 -sb ln:-14,1 -sl u:30,80 -st u:100000,10000000 -si u:1,6 -so ln:0,0.1 -gb u:-25,-15 -gt f:gb -ld ln:gb,0.4 -lb f:ld -lt ln:gt,0.4 -lk 1 -sp ln:12,0.5 -sg f:1 -hs ln:1.5,1 -hl ln:1.2,1 -hg ln:1.4,1 -su e:10000000 -o SP_Tree_generate -om 1 -v 1 -od 1 -oc 1 -op 1 -ol 1 -cs $RANDOM 
	cp SP_Tree_generate/1/s_tree.trees  s_tree.trees 

	python simphy_wgd.py -spfile=s_tree.trees -auxfile=myaux.txt -sout=s_tree_modded.trees
	python simphy_wgd.py -spfile=s_tree_modded.trees -auxfile=myaux2.txt -sout=s_tree_modded2.trees
	s_tree_modded=$(cat s_tree_modded2.trees)
	./simphy -sb f:0.000001 -gt f:0.0 -gg f:0.0 -gp f:0.0 -ld f:0.000000000000000000001 -lb f:0.000000000000000000001 -lt f:0.0 -lg f:0.0 -rs 1 -rl f:99 -rg 1 -o $dir_name -sp f:2 -su f:0.00001 -sg f:1 -sl U:50,110 -st f:1000000 -om 1 -v 3 -od 1 -op 1 -oc 1 -on 1 -ol 1 -cs $RANDOM -S $s_tree_modded
	
	cd $dir_name
	newname="sim_${i}"
	rm -r $newname
	mv 1 $newname
	cd ..

	for (( j=1; j<=99; j++ ))
	do
    		ii=$(printf "%02d" $j)
		filename="${ii}.mapsl"
		filename2="${ii}l1g.maplg"
		filename3="g_trees${ii}.trees"
		echo $filename
		python simphy_wgd.py -mode=remap -auxfile=myaux2.txt -slfile=$dir_name/$newname/$filename -lgfile=$dir_name/$newname/$filename2 -genetreefile=$dir_name/$newname/$filename3 > out.txt
	done
	
	for (( j=1; j<=99; j++ ))
	do
    		ii=$(printf "%02d" $j)
		filename="${ii}.mapsl.modded"
		filename2="${ii}l1g.maplg.modded"
		filename3="g_trees${ii}.trees.modded"
		echo $filename
		python simphy_wgd.py -mode=remap -auxfile=myaux.txt -slfile=$dir_name/$newname/$filename -lgfile=$dir_name/$newname/$filename2 -genetreefile=$dir_name/$newname/$filename3 > out.txt
	done


	cd $dir_name
	cp ${dir_name}.command $newname/${newname}.command
	cd $newname 
	cat g_trees*.modded.modded >> all_genetrees.txt
	echo "gene trees are combined!"

	cd ..
	cd ..
	cp myaux.txt $dir_name/$newname/myaux.txt
	cp myaux2.txt $dir_name/$newname/myaux2.txt
	cp s_tree.trees $dir_name/$newname/source_s_tree.trees
	cp s_tree_modded.trees $dir_name/$newname/s_tree_modded.trees
	cp s_tree_modded2.trees $dir_name/$newname/s_tree_modded2.trees

	
	#maxheight=$(head -n 1 "$dir_name/$newname/maxheight.txt")
	maxheight="-"
	
	if [ -s $dir_name/$newname/s_tree.trees ]; then

		python post-order-labeling.py $dir_name/$newname/source_s_tree.trees $dir_name/$newname/s_tree.newick
		python map_gene_trees_twoWGD.py $dir_name/$newname/all_genetrees.txt $dir_name/$newname/all_genetrees_edited.txt $dir_name/$newname

	fi
done


#apply losses

lossrates=("decNum" "decBoth" "fix")

for (( i=1; i<=$Numbersimulation; i++ ))
do	
	newname="sim_${i}"
	if [ -s $dir_name/$newname/s_tree.trees ]; then
		for loss in "${lossrates[@]}"; do
			if [ $loss == "decNum" ]; then
				out_name="applied_loss_decNum_all_genetrees_edited.txt"
				python apply_losses_on_simphy.py $dir_name/$newname/$out_name $dir_name/$newname/all_genetrees_edited.txt 1 0 0 1
			elif [ $loss == "decBoth" ]; then
				out_name="applied_loss_decBoth_all_genetrees_edited.txt"
				python apply_losses_on_simphy.py $dir_name/$newname/$out_name $dir_name/$newname/all_genetrees_edited.txt 1 1 0 0
			elif [ $loss == "fix" ]; then
				out_name="applied_loss_fix_all_genetrees_edited.txt"
				python apply_losses_on_simphy.py $dir_name/$newname/$out_name $dir_name/$newname/all_genetrees_edited.txt 1 0 1 1

			fi
		done

	else

		stats=-1
	
	fi
	echo "$newname done!"
done
