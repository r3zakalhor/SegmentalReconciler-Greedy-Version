#!/bin/sh
declare -i bool
declare -i nbgenetrees
declare -a duparr
nbgenetrees=1000
bool=1
duparr=(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0) # size of 120
for (( i=1; i<=$nbgenetrees; i++ ))
do
#echo ${#i}
j=$i
if ((j < 10)); then
    j="00$j"
    #echo $j
elif ((j < 100)); then
    j="0$j"
fi
while read line; do
    for word in $line; do
        
        if (( $bool == 0 )); then
		#echo $word
		int=`echo $word | cut -d "'" -f 2`
  		#echo $int
		let duparr[$int]++
		#echo ${duparr[$int]}
		bool=1
	
	fi

	if [ $word = "Dup" ]
	then
		#echo "Yes!"
		bool=0
	
	fi
    done
done <"$j.mapsl"
done
maxheight=0
for i in ${!duparr[@]}; do
  if [ $maxheight -lt ${duparr[$i]} ]
  then
      maxheight=${duparr[$i]}
      echo "$maxheight"
  fi
  echo "nb of dup for species $i is ${duparr[$i]}" >> nb_dups_species.txt
done

echo "$maxheight" >> maxheight.txt