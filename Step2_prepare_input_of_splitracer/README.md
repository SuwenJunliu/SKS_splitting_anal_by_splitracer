# Step2. prepare_input_of_splitracer

The input of the splitracer consists of the following components

- station.lst file
- data_request_$Station_Name$_$Network_Name$ file
- channel_info.txt file

All those file are text files, and the demo files are in ./demo_input

## station.lst file

The file which contains the station location information with the following format

Station_name Network_name latitude longitude


## data_request file

The file which contains the SAC files and earthquake events information with the following format

year month day hour minute second evt_latitude evt_longitude gcarc baz depth mag loc filepath

NOTE: the filepath is the `full path` without the prefix, i.e.,  D:\NECESSArrayEvts\20090924071622\20090924071622.YP.NE11.BHE means that there are three sacfiles 20090924071622.YP.NE11.BH[N,E,Z].sac

## channel_info.txt file

The file which contains the station 3 component information with the following format

Network | Station | Location | Channel | Latitude | Longitude | Elevation | Depth | Azimuth | Dip | SensorDescription | Scale | ScaleFreq | ScaleUnits | SampleRate | Starttime | EndTime
YP NE11 - BHE 42.85992 124.05263  0 90 0 NoDescrip 1 1 Count 100 1979-01-01T00:00:00 2099-01-01T00:00:00 
YP NE11 - BHN 42.85992 124.05263  0 0 0 NoDescrip 1 1 Count 100 1979-01-01T00:00:00 2099-01-01T00:00:00 
YP NE11 - BHZ 42.85992 124.05263  0 0 -90 NoDescrip 1 1 Count 100 1979-01-01T00:00:00 2099-01-01T00:00:00 

NOTE: include header

## Script
There are serval script to generate this three files