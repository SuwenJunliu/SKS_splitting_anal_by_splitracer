#!/bin/bash


## para
minMag=5.5
minGcarc=85
maxGcarc=125
##


inDir="/mnt/d/NECESSArrayEvts/"
outDir="/mnt/d/NECESSArrayEvts/output"


cd $inDir
ls -d 20* | sort -g  > eventTemp
while read eventFolder
do
	echo "deal with station folder $eventFolder"
	cd $eventFolder
	ls *.sac  | awk -F"." '{print $3}' | sort -g | uniq > stationNameTemp
	

	while read stationName
	do
	
		# get the info by saclst
		echo "deal with $stationName"
		mag=`saclst mag f *$stationName* | awk '{print $2}' | head -n 1`
		gcarc=`saclst gcarc f *$stationName* | awk '{print $2}' | head -n 1`	

		if [[ `echo "$mag < $minMag" | bc` -eq 1 ]];then
			echo "mag $mag is less than $minMag, skip"
			break
		fi
		if [[ `echo "$gcarc < $minGcarc" | bc` -eq 1 ]] || [[ `echo "$gcarc > $maxGcarc" | bc` -eq 1 ]];then
			echo "gcarc $gcarc is less or greater than $minGcarc or $maxGcarc, skip"
			break
		fi
		
		filePrefix=`ls *$stationName* | head -n 1 | awk -F"." '{print $1}'`
		network=`saclst knetwk f *$stationName* | awk '{print $2}' | head -n 1`
		stla=`saclst stla f *$stationName* | awk '{print $2}' | head -n 1`
		stlo=`saclst stlo f *$stationName* | awk '{print $2}' | head -n 1`
		ss=`saclst kztime f *$stationName* | awk '{print $2}' | head -n 1 | awk -F":" '{print $3}'`

		#echo "prefix $filePrefix"
		yyyy=`echo ${filePrefix:0:4}`
		mm=`echo ${filePrefix:4:2}`
		dd=`echo ${filePrefix:6:2}`
		hh=`echo ${filePrefix:8:2}`
		mmin=`echo ${filePrefix:10:2}`
		#ss=`echo ${filePrefix:12:2}`
		
		evla=`saclst evla f *$stationName* | awk '{print $2}' | head -n 1`
		evlo=`saclst evlo f *$stationName* | awk '{print $2}' | head -n 1`
		baz=`saclst baz f *$stationName* | awk '{print $2}' | head -n 1`
		depth=`saclst evdp f *$stationName* | awk '{print $2}' | head -n 1`

		echo "$yyyy $mm $dd $hh $mmin $ss $evla $evlo $gcarc $baz $depth $mag - D:\\NECESSArrayEvts\\$eventFolder\\$eventFolder.${network}.${stationName}" >> $outDir/data_request_${stationName}_$network	
		echo "good station. $stationName"
	
	done < stationNameTemp
	
	rm stationNameTemp

	cd ..
	pwd
done < eventTemp

rm eventTemp


