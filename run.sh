#!/bin/bash

#	char *params_name = argv[1];
#       char *matrix_weather_name = argv[2];
#       char *days_harvest_name = argv[3];
#       int ndays = atoi(argv[4]);
#       int nout = atoi(argv[5]);
#       int y = atoi(argv[6]);

params_name="params.csv"
matrix_weather_name="weather.csv"
days_harvest_name="days_harvest.csv"
ndays=489
year_start=1999
doy_start=227

./basgra $params_name $matrix_weather_name $days_harvest_name $ndays $doy_start $year_start
