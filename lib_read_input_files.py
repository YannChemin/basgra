import numpy as np

# ModVege is using 4 parts:
# Green Vegetative  (GV)
# Green Reproductive(GR)
# Dead Vegetative   (DV)
# Dead Reproductive (DR)

#########################################################
# Definition of input parameters
#########################################################
# station name              saerheim_summer_00_gri
# log10clvi                 1
# log10cresi                0
# log10crti                 1
# csti                      0
# log10laii                 0.1
# pheni                     0.01
# tiltoti                   800
# frtilgi                   0
# lt50i                     -4.7894
# claiv                     0.5
# cocresmx                  0.141155
# cstavm                    0.229575
# daylb                     0.539664
# daylp                     0.632348
# dlmxge                    0.936974
# fslamin                   0.686948
# fsmax                     0.693
# hagere                    0.8
# k                         0.607039
# laicr                     3.79569
# laieft                    0.2
# laitil                    0.566455
# lfwidg                    0.0085207
# lfwidv                    0.00491875
# nellvm                    2.09178
# phencr                    0.495209
# phy                       63.0532
# rdrsco                    0.0712247
# rdrsmx                    0.06
# rdrtem                    0.00102573
# rgenmx                    0.0108797
# rootdm                    0.6
# rrdmax                    0.012
# rubisc                    5.7803
# shape                     0.538907
# simax1t                   0.00450449
# slamax                    0.06
# tbase                     3.61083
# tcres                     1.88986
# toptge                    12.6178
# tranco                    8
# yg                        0.841767
# lat                       58.46
# wci                       0.3
# fwcad                     0.011363636
# fwcwp                     0.056818182
# fwcfc                     0.681818182
# fwcwet                    1
# wcst                      0.44
# wpoolmax                  50
# dparam                    0.0031796
# fgas                      0.4
# fo2mx                     0.21
# gamma                     65
# hparam                    0.0055925
# krdranaer                 0.2
# kresphard                 0.1
# krsr3h                    1
# krtotaer                  2
# ksnow                     0.035
# lambdasoil                172800
# ldt50a                    1.3403
# ldt50b                    -2.1128
# lt50mn                    -26.6839
# lt50mx                    -4.7894
# ratedmx                   2
# rehardredday              145
# rhonewsnow                100
# rhopack                   0.02
# swret                     0.1
# swrf                      0.01
# thardmx                   14.7052
# tmeltfreeze               0
# trainsnow                 0.01
# tsurfdiff                 0.62279
# kluetilg                  0.5
# frtilgg1i                 0.1
# daylg1g2                  0.6
# rgrtg1g2                  0.9
# rdrtmin                   0.01
# tvern                     20


def read_params(input_params_csv):
    arr = np.genfromtxt(input_params_csv, delimiter=",", dtype=float, skip_header=1, usecols=(-1))
    return(arr)

#Test only
#input_params_csv='params.csv'
#print(read_params(input_params_csv))

#########################################################
# Definition of input parameters
#########################################################
# year          Year
# doy	        Day Of Year
# rdd           Ground Radiation (MJ.m-2)
# tmin          Minimum temperature (degree celsius)
# tmax          Maximum temperature (degree celsius)
# vp            Vapour Pressure (kPa) 
# wn            Wind Speed (m/s)
# rain          Precipitation (mm/day)
# eact          Actual Evaporation from Remote Sensing (mm/day)
# tact          Actual Transpiration from Remote Sensing (mm/day)
# lai           Leaf Area Index from Remote Sensing (cm2/cm2) 
# gcut          Grass cut event from Remote Sensing temporal Coherence (0/1)


def read_weather(input_weather_csv):
    """
    @params input_weather_csv the input file named weather.csv
    @params arr the returning array of [year, doy, rdd, tmin, tmax, vp, wn, rain, eact, tact, lai, gcut]
    """
    arr = np.genfromtxt(input_weather_csv, delimiter=",", dtype=float, skip_header=1)
    return(arr)


def read_calval(input_calval_csv):
    """
    @params input_calval_csv the input file named calval.csv
    @params arr the returning array of [YEAR,DOY,RES,DM,LAI,LERG,NELLVG,RLEAF,SLA,TILTOT,FRTILG,CP]

    """
    arr = np.genfromtxt(input_calval_csv, delimiter=",", dtype=float, skip_header=1)
    return(arr)


#Test only
#input_weather_csv='weather.csv'
#print(read_weather(input_weather_csv))

#Test only
#input_calval_csv='calval.csv'
#print(read_calval(input_calval_csv))
