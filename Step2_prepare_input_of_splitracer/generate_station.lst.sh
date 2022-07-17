#!/bin/bash

outFile="./station.lst"
while read line
do
	echo $line | awk '{print $1,"YP",$2,$3}' >> $outFile
done < YP.coords
