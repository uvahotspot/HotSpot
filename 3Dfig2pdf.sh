#!/bin/bash

ary=()

while read data; do
    ary=( "${ary[@]}" "$data" )
done

COUNT="0"

temp=()
echo -e "\nLoading Files:"
for i in ${ary[*]}; do
    echo "$(($COUNT+1)): $i"
    let COUNT++;
    name=${i%.fig}
    rand="$name".$RANDOM.$$.ps
    temp=( "${temp[@]}" "$rand")
    fig2dev -L ps "$i" 1>"$rand"
done

COUNT="0"
echo -e "\nCreated Temp Files:"
for j in ${temp[*]}; do
    echo "$(($COUNT+1)): $j"
    let COUNT++;
done

merge=mergedps.$RANDOM.$$.ps
psmerge -o$merge  ${temp[*]}
echo "$(($COUNT+1)): $merge"

ps2pdf $merge $1
echo -e "\nCreated PDF: $1"



COUNT="0"
echo -e "\nRemoving Temp Files:"
for k in ${temp[*]}; do
    echo "$(($COUNT+1)): $k"
    let COUNT++;
    rm $k
done

echo "$(($COUNT+1)): $merge"
rm $merge

let COUNT=0;
echo -e "Removing .Fig Files:"
deleted=()
for i in ${ary[*]}; do
    delete="1"
    for name in ${deleted[*]}; do
        if [ "$i" == "$name" ];
        then
            let delete="0"
        fi
    done
    if [ "$delete" == "1" ];then
        rm $i
        deleted+=($i)
        echo "$(($COUNT+1)): '$i'"
    fi
    let COUNT++;

done

