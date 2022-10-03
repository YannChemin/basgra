# Read input parameters
from lib_read_input_files import *
# Internal basgra functions
from lib_basgra import *
# basgra model as a function
from basgra import *
# Read arguments given to this script
import argparse
parser = argparse.ArgumentParser(description="Runs the BASGRA simulation")
parser.add_argument("params_fname", default='params.csv', help="The CSV file holding the input parameters")
parser.add_argument("weather_fname", default='weather.csv', help="The CSV file holding the weather data, the RS data (Evap, tran, LAI & cut dates)")
parser.add_argument("startdoy", default=1, type=int, help="The DOY at which the simulation starts")
parser.add_argument("enddoy", default=365, type=int, help="The DOY at which the simulation stops")
parser.add_argument("startyear", default=2002, type=int, help="The YEAR at which the simulation starts")
parser.add_argument("endyear", default=2002, type=int, help="The YEAR at which the simulation stops")
args = parser.parse_args()

# Input definition 
#params_fname = 'params.csv'
#weather_fname = 'weather.csv'
#startdoy = 1
#enddoy = 365
#startyear = 2002
#endyear = 2002

# Create out_dates list
import datetime
out_dates = []
d = datetime.date(args.startyear,1,1)+datetime.timedelta(days=args.startdoy)
while d <= datetime.date(args.endyear,1,1)+datetime.timedelta(days=args.enddoy):
    out_dates.append(d)
    d += datetime.timedelta(days=1)


# Read Params file (params.csv) and create array
params = read_params(args.params_fname)
# Read Weather file (weather.csv) and create array
weather = read_weather(args.weather_fname)

def run_basgra(params, weather, startdoy, enddoy, startyear, endyear):
    # Run BASGRA
    #ph, evp, trn, tl_v, tl_n, tl_h, gv_b, dv_b, h_b = basgra(params, weather, startdoy, enddoy, startyear, endyear)
    nell_vg, laim, ph, evp, trn, tl_v, tl_n, tl_h, gv_b, dv_b, h_b = basgra(params, weather, startdoy, enddoy, startyear, endyear)

    # Input CALVAL #### ONLY FOR DEBUG !
    calval_fname = 'calval.csv'
    # Read CALVAL #### ONLY FOR DEBUG !
    calval = read_calval(calval_fname)
    #YEAR,DOY,RES,DM,LAI,LERG,NELLVG,RLEAF,SLA,TILTOT,FRTILG,CP
    # calendar & weather
    year    = calval[:,0]  # YYYY
    doy     = calval[:,1]  # DOY
    res     = calval[:,2]  # RES ?
    dm      = calval[:,3]  # DM [kg/ha?]
    lai     = calval[:,4]  # LAI [-]
    lerg    = calval[:,5]  # LERG: Elongation rate of leaves on elongating tillers (m/day)
    nellvg  = calval[:,6]  # NELLVG: Number of growing leaves per elongating tiller
    rleaf   = calval[:,7]  # RLEAF: reductions in leaf appearance rate
    sla     = calval[:,8]  # SLA: specific leaf area 
    tiltot  = calval[:,9]  # Tillers Total 
    frtilg  = calval[:,10] # FRTILG ? (a single value of 0.8 in a year...)
    cp      = calval[:,11] # CP ?

    #create the Harvesting list
    harvested = [0] * len(dm)
    harvested[150] = dm[150]-dm[151]
    for i in range(151,len(harvested),1):
        harvested[i] = harvested[150]


    #resize out_dates to actual output size
    out_d = out_dates[:len(gv_b)]
    import numpy as np
    import matplotlib.pyplot as plt

    phenrf = np.multiply(np.subtract(1,ph),1.981017)
    phenrf_re = np.divide(nell_vg,2.09178)

    plt.figure(figsize=(15,7))

    plt.subplot(331)
    plt.plot(out_d,np.multiply(10,gv_b),'m-',label="gv_b")
    plt.plot(out_d,dm[:len(out_d)],'g-',label="gv_b ref")
    plt.title('Green Vegetative biomass (kg DM/ha)')
    plt.legend()
    plt.grid()

    plt.subplot(332)
    plt.plot(out_d,np.add(evp,trn),'c-',label="ET")
    plt.plot(out_d,evp,'b-',label="evaporation")
    plt.plot(out_d,trn,'g-',label="transpiration")
    plt.title('Evaporation and Transpiration (mm/d)')
    plt.legend()
    plt.grid()

    plt.subplot(333)
    plt.plot(out_d,np.add(tl_v,tl_h),'r-',label="tillers")
    plt.plot(out_d,tiltot[:len(out_d)],'b-',label="tillers ref")
    plt.title('Tillers')
    plt.legend()
    plt.grid()

    plt.subplot(334)
    plt.plot(out_d,dv_b,'g-',label="dv_b")
    #plt.plot(out_d,out_dvb,'b-',label="out_dvb")
    plt.title('Dead Vegetative biomass (kg DM/ha)')
    plt.legend()
    plt.grid()

    plt.subplot(335)
    plt.plot(out_d,ph,'b-',label="phen")
    plt.plot(out_d,phenrf,'g-',label="phenrf ext. calc.")
    plt.plot(out_d,phenrf_re,'m-',label="phenrf rev Engg")
    plt.title('Phenology')
    plt.legend()
    plt.grid()

    plt.subplot(336)
    plt.plot(out_d,tl_n,'b-',label="tilg1: Non-elong generative")
    plt.plot(out_d,tl_h,'g-',label="tilg2: Elong. generative")
    plt.title('Tillers Count (tillers/m2)')
    plt.legend()
    plt.grid()

    plt.subplot(337)
    plt.plot(out_d,laim,'m-',label="LAI")
    plt.plot(out_d,lai[:len(out_d)],'g-',label="LAI ref")
    plt.title('LAI')
    plt.legend()
    plt.grid()

    # Harvested Biomass Plot
    plt.subplot(338)
    plt.plot(out_d,h_b,'m-',label="h_b")
    plt.plot(out_d,harvested[:len(out_d)],'g-',label="h_b ref")
    plt.title('Harvested biomass (kg DM/ha)')
    plt.legend()
    plt.grid()

    plt.subplot(339)
    plt.plot(out_d,nell_vg,'m-',label="nellvg")
    plt.plot(out_d,nellvg[:len(out_d)],'g-',label="nellvg ref")
    #plt.plot(out_d,tl_h,'c-',label="tillers harvestable")
    plt.title('NELL VG')
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.show()

# run the main function
run_basgra(params, weather, args.startdoy, args.enddoy, args.startyear, args.endyear)
