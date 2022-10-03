#!/bin/bash

# Input definition 
params_fname='params.csv'
weather_fname='weather.csv'
startdoy=1
enddoy=365
startyear=2002
endyear=2002

python run_basgra.py $params_fname $weather_fname $startdoy $enddoy $startyear $endyear

