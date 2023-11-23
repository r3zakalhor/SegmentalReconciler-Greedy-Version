#!/bin/sh
declare -i bool
declare -i nbgenetrees
declare -a duparr
nbgenetrees=100
bool=1

Numbersimulation=100

dir_name="s4_100"
cd $dir_name
mkdir nb_dup_species  # Create the directory if it doesn't exist

str="species id,nb dup"




write_csv(){
    echo \"$1\",\"$2\" >> $stats_out
}

for (( i=1; i<=$Numbersimulation; i++ )); do

    
    
    name="nb_dup_sim_${i}.csv"
    stats_out=nb_dup_species/$name
    echo $str > $stats_out

    newname="sim_${i}"
    cd $newname 

    duparr=(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0) # size of 120

    for (( j=1; j<=$nbgenetrees; j++ )); do
        if ((j < 10)); then
            padded_j="0$j"
        elif ((j < 100)); then
            padded_j="$j"
        fi

        while read line; do
            for word in $line; do
                if (( $bool == 0 )); then
                    int=`echo $word | cut -d "'" -f 2`
                    let duparr[$int]++
                    bool=1
                fi

                if [ $word = "Dup" ]; then
                    bool=0
                fi
            done
        done <"$padded_j.mapsl"
    done
    cd ..
	
    echo "$i simulation counted!"    
	
    for k in ${!duparr[@]}; do
        write_csv $k ${duparr[$k]}
    done

    

done
