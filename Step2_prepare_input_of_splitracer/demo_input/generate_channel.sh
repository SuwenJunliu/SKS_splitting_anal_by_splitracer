#!/bin/bash

#outFile="./channel_info.txt"
output_path="/mnt/d/NECESSArrayEvts/output"
#echo "Network | Station | Location | Channel | Latitude | Longitude | Elevation | Depth | Azimuth | Dip | SensorDescription | Scale | ScaleFreq | ScaleUnits | SampleRate | Starttime | EndTime"
while read line
do
	station_name=`echo $line | awk '{print $1}'`
	network_name="YP"
	mkdir ${output_path}/${station_name}_${network_name}
	echo "Network | Station | Location | Channel | Latitude | Longitude | Elevation | Depth | Azimuth | Dip | SensorDescription | Scale | ScaleFreq | ScaleUnits | SampleRate | Starttime | EndTime" >> $output_path/${station_name}_${network_name}/channel_info.txt
	echo $line | awk '{print "YP",$1,"-","BHE",$2,$3,$5,"0","90","0","NoDescrip","1","1","Count","100","1979-01-01T00:00:00","2099-01-01T00:00:00 "}' >> $output_path/${station_name}_${network_name}/channel_info.txt
	echo $line | awk '{print "YP",$1,"-","BHN",$2,$3,$5,"0","0","0","NoDescrip","1","1","Count","100","1979-01-01T00:00:00","2099-01-01T00:00:00 "}' >> $output_path/${station_name}_${network_name}/channel_info.txt
	echo $line | awk '{print "YP",$1,"-","BHZ",$2,$3,$5,"0","0","-90","NoDescrip","1","1","Count","100","1979-01-01T00:00:00","2099-01-01T00:00:00 "}' >> $output_path/${station_name}_${network_name}/channel_info.txt

done < YP.coords
