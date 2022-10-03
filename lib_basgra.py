
from math import sin,cos,tan,asin,atan,pi,sqrt,exp,pi

####################################################################################################
# module environment
####################################################################################################

def rainsnowsurfacepool(davtmp, trainsnow, rain, bias, ampl, freq, doy, tmeltfreeze, melt, swret, swrf, drystor, snowmelt, wetstor, refreeze, pwater, wavail, sdepth, delt, rhopack, lai, fdepth, poolinfillimit, pinfil, wpoolmax, wapl, waps, runon, poolvolremain, poolinfil, frate, tsurf, lambdaice, rhowater, latentheat, eta, freezepl, infil, packmelt, pooldrain, thawps, wremain):
    """
    Rain Snow Surface Pool

    @param davtmp Daily average temperature (Celsius)
    @param trainsnow Temperature below which precipitation is snow (Celsius)
    @param rain (mm/day)
    @param bias Average snow melting rate at 1 degree above 'TmeltFreeze' (mm/degreeCelsius/day)
    @param ampl Intra-annual amplitude snow melt at 1 degree > 'TmeltFreeze' (mm/degreeCelsius/day)
    @param freq Unit conversion time to annual cycle (radians/day)
    @param doy Day Of Year (1-366)
    @param tmeltfreeze Temperature above which snow melts (Celsius)
    @param melt Potential snow melt rate per degree above TmeltFreeze (mm/C/day)
    @param swret Liquid water storage capacity of snow (mm/mm/day)
    @param swrf Maximum refreezing rate per degree below 'TmeltFreeze' (mm/day/Celsius)
    @param drystor Snow amount as SWE (Soil Water Equivalent)
    @param snowmelt Snow melting (mm/day)
    @param wetstor Liquid water in snow (mm)
    @param refreeze Freezing of liquid water stored in snow (mm/day)
    @param pwater Potential water available from rain (mm/day)
    @param wavail Liquid water from rain, snow melt and storage in snow (mm/day)
    @param sdepth Snow depth (m)
    @param delt Model time step (day)
    @param rhopack Relative packing rate of snow (-/day)
    @param lai Leaf Area Index (m2 leaf / m2)
    @param fdepth Soil frost layer depth (m)
    @param poolinfillimit Soil frost depth limit for water infiltration (m)
    @param pinfil Infil of Pool (Wsupply - RNINTC) (mm/day)
    @param wpoolmax Maximum pool water (liquid plus ice) (mm)
    @param wapl Pool water amount: liquid (mm)
    @param waps Pool water amount: solid (=ice) (mm)
    @param runon Water in excess of what can infiltrate the soil (mm/day)
    @param poolvolremain Unused capacity of surface pool (mm/day)
    @param poolinfil Water flow to pool from other sources than ice thawing (mm/day)    
    @param frate Rate of increase of frost layer depth (m/day)
    @param tsurf Soil surface temperature (C)
    @param lambdaice Thermal conductivity of ice (J/m/K/day)
    @param rhowater Density of water (kg/m3)
    @param latentheat Latent heat of water fusion (J/kg)
    @param eta Compound parameter (= LAMBDAice/(RHOwater*LatentHeat)) (m2/K/day)
    @param freezepl Freezing rate of pool water (mm/day)
    @param infil Water flow into soil from precipitation and snow melt (mm/day)
    @param packmelt Loss of snow height by packing and by melting (m/day)
    @param pooldrain Water flow from pool to soil (mm/day)
    @param thawps Rate of surface ice thawing (mm/day)
    @param wremain Liquid water staying in snow pack (mm/day)
    @result pwater Potential water available from rain (mm/day)
    @result psnow Potential snow available from rain (mm/day)
    @result snowmelt Snow melting (mm/day)
    @result wmaxstore Liquid water storage capacity of the snowpack (mm/day)
    @result refreeze Freezing of liquid water stored in snow (mm/day)
    @result staywet Liquid water in snow remaining liquid (mm/day)
    @result wavail Liquid water from rain, snow melt and storage in snow (mm/day)
    @result wremain Liquid water staying in snow pack (mm/day)
    @result wsupply Liquid water not staying in snow pack (mm/day)
    @result density Density of snow (m3/kg)
    @result packmelt Loss of snow height by packing and by melting (m/day)
    @result rnintc Interception of precipitation by the canopy (mm/day)
    @result pinfil Infil of Pool (Wsupply - RNINTC) (mm/day)
    @result runon Water in excess of what can infiltrate the soil (mm/day)
    @result freezepl Freezing rate of pool water (mm/day)
    @result thawps Rate of surface ice thawing (mm/day)
    @result pooldrain Water flow from pool to soil (mm/day)
    """
    pwater, psnow = precform(davtmp, trainsnow, rain)
    snowmelt, wmaxstore, refreeze, staywet, wavail, wremain, wsupply, density, packmelt = watersnow(bias, ampl, freq, doy, davtmp, tmeltfreeze, drystor, delt, melt, swret, swrf, wetstor, refreeze, pwater, wavail, sdepth, rhopack)
    rnintc = min( wsupply, 0.25*lai)
    pinfil = wsupply - rnintc
    runon = infilrunon(fdepth, poolinfillimit, pinfil)
    freezepl, thawps, pooldrain = surfacepool(wpoolmax, wapl, waps, runon, poolvolremain, poolinfil, delt, fdepth, poolinfillimit, frate, tsurf, lambdaice, rhowater, latentheat, eta)
    return(pwater, psnow, snowmelt, wmaxstore, refreeze, staywet, wavail, wremain, wsupply, density, packmelt, rnintc, pinfil, runon, freezepl, thawps, pooldrain)

def precform(davtmp, trainsnow, rain):
    """
    Partition rainfall volume into either potential water or potential snow depending on temperature

    @param davtmp Daily average temperature (Celsius)
    @param trainsnow Temperature below which precipitation is snow (Celsius)
    @param rain (mm/day)
    @result pwater Potential water available from rain (mm/day)
    @result psnow Potential snow available from rain (mm/day)
    """
    if (davtmp > trainsnow):
        pwater = rain
        psnow  = 0.
    else:
        pwater = 0.
        psnow  = rain

    return(pwater, psnow)

def watersnow(bias, ampl, freq, doy, davtmp, tmeltfreeze, drystor, delt, melt, swret, swrf, wetstor, refreeze, pwater, wavail, sdepth, rhopack):
    """
    Water Snow

    @param bias Average snow melting rate at 1 degree above 'TmeltFreeze' (mm/degreeCelsius/day)
    @param ampl Intra-annual amplitude snow melt at 1 degree > 'TmeltFreeze' (mm/degreeCelsius/day)
    @param freq Unit conversion time to annual cycle (radians/day)
    @param doy Day Of Year (1-366)
    @param davtmp Daily average temperature (Celsius)
    @param tmeltfreeze Temperature above which snow melts (Celsius)
    @param drystor Snow amount as SWE (Soil Water Equivalent)
    @param delt Model time step (day)
    @param melt Potential snow melt rate per degree above TmeltFreeze (mm/C/day)
    @param swret Liquid water storage capacity of snow (mm/mm/day)
    @param swrf Maximum refreezing rate per degree below 'TmeltFreeze' (mm/day/Celsius)
    @param wetstor Liquid water in snow (mm)
    @param refreeze Freezing of liquid water stored in snow (mm/day)
    @param pwater Potential water available from rain (mm/day)
    @param wavail Liquid water from rain, snow melt and storage in snow (mm/day)
    @param sdepth Snow depth (m)
    @param rhopack Relative packing rate of snow (-/day)
    @result snowmelt Snow melting (mm/day)
    @result wmaxstore Liquid water storage capacity of the snowpack (mm/day)
    @result refreeze Freezing of liquid water stored in snow (mm/day)
    @result staywet Liquid water in snow remaining liquid (mm/day)
    @result wavail Liquid water from rain, snow melt and storage in snow (mm/day)
    @result wremain Liquid water staying in snow pack (mm/day)
    @result wsupply Liquid water not staying in snow pack (mm/day)
    @result density Density of snow (m3/kg)
    @result packmelt Loss of snow height by packing and by melting (m/day)
    """
    snowmelt, wmaxstore         = snowmeltwmaxstore(bias, ampl, freq, doy, davtmp, tmeltfreeze, drystor, delt, melt, swret)
    refreeze, staywet           = wetstordynamics(swrf, tmeltfreeze, davtmp, wetstor, delt)
    wavail, wremain, wsupply    = liquidwaterdistribution(staywet, snowmelt, pwater, wavail, wmaxstore)
    density                     = snowdensity(drystor, wetstor, sdepth)
    packmelt                    = snowdepthdecrease(sdepth, delt, rhopack, snowmelt, density)
    return(snowmelt, wmaxstore, refreeze, staywet, wavail, wremain, wsupply, density, packmelt)

def snowmeltwmaxstore(bias, ampl, freq, doy, davtmp, tmeltfreeze, drystor, delt, melt, swret):
    """
    Snow Melt Water Max Storage

    @param bias Average snow melting rate at 1 degree above 'TmeltFreeze' (mm/degreeCelsius/day)
    @param ampl Intra-annual amplitude snow melt at 1 degree > 'TmeltFreeze' (mm/degreeCelsius/day)
    @param freq Unit conversion time to annual cycle (radians/day)
    @param doy Day Of Year (1-366)
    @param davtmp Daily average temperature (Celsius)
    @param tmeltfreeze Temperature above which snow melts (Celsius)
    @param drystor Snow amount as SWE (Soil Water Equivalent)
    @param delt Model time step (day)
    @param melt Potential snow melt rate per degree above TmeltFreeze (mm/C/day)
    @param swret Liquid water storage capacity of snow (mm/mm/day)
    @result snowmelt Snow melting (mm/day)
    @result wmaxstore Liquid water storage capacity of the snowpack (mm/day)
    """
    melt = bias + ampl * sin( freq * (doy-(174.-91.)) )
    if (davtmp > tmeltfreeze):
        snowmelt = max( 0., min( drystor/delt, melt*(davtmp-tmeltfreeze) ))
    else:
        snowmelt = 0.

    wmaxstore = drystor * swret
    return(snowmelt, wmaxstore)

def wetstordynamics(swrf, tmeltfreeze, davtmp, wetstor, delt):
    """
    Wet Storage Dynamics

    @param swrf Maximum refreezing rate per degree below 'TmeltFreeze' (mm/day/Celsius)
    @param tmeltfreeze Temperature above which snow melts (Celsius)
    @param davtmp Daily average temperature (Celsius)
    @param wetstor Liquid water in snow (mm)
    @param delt Model time step (day)
    @result refreeze Freezing of liquid water stored in snow (mm/day)
    @result staywet Liquid water in snow remaining liquid (mm/day)
    """
    refreezemax = swrf * (tmeltfreeze-davtmp)
    if ((wetstor > 0) and (davtmp < tmeltfreeze)):
        refreeze = min(wetstor/delt, refreezemax)
    else:
        refreeze = 0.

    staywet = wetstor/delt - refreeze
    return(refreeze, staywet)

def liquidwaterdistribution(staywet, snowmelt, pwater, wavail, wmaxstore):
    """
    Liquid Water Distribution

    @param staywet Liquid water in snow remaining liquid (mm/day)
    @param snowmelt Snow melting (mm/day)
    @param pwater Potential water vailable from rain (mm/day)
    @param wavail Liquid water from rain, snow melt and storage in snow (mm/day)
    @param wmaxstore Liquid water storage capacity of the snowpack (mm/day)
    @result wavail Liquid water from rain, snow melt and storage in snow (mm/day)
    @result wremain Liquid water staying in snow pack (mm/day)
    @result wsupply Liquid water not staying in snow pack (mm/day)
    """
    wavail  = staywet + snowmelt + pwater
    wremain = min(wavail,wmaxstore)
    wsupply = wavail - wremain
    return(wavail, wremain, wsupply)

def snowdensity(drystor, wetstor, sdepth):
    """
    Snow Density

    @param drystor Snow amount as SWE (Soil Water Equivalent)
    @param wetstor Liquid water in snow (mm)
    @param sdepth Snow depth (m)
    @result density Snow Density (m3/kg)
    """
    swe = drystor + wetstor
    if (sdepth > 0.):
        return min(480., swe/sdepth)
    else :
        return 0.

def snowdepthdecrease( sdepth, delt, rhopack, snowmelt, density):
    """
    Snow Depth Decrease

    @param sdepth Snow depth (m)
    @param delt Model time step (day)
    @param rhopack Relative packing rate of snow (-/day)
    @param snowmelt Snow melting (mm/day)
    @param density Density of snow (m3/kg)
    @result snowdepth Decrease in snow depth
    """
    if (sdepth > 0. and density > 0):
        return max(0.,min( sdepth/delt, sdepth*rhopack - snowmelt/density ))
    else :
        return 0.

def infilrunon(fdepth, poolinfillimit, pinfil):
    """
    Infiltration Run-On

    @param fdepth Soil frost layer depth (m)
    @param poolinfillimit Soil frost depth limit for water infiltration (m)
    @param pinfil Infil of Pool (Wsupply - RNINTC) (mm/day)
    @result runon Water in excess of what can infiltrate the soil (mm/day)
    """
    if (fdepth <= poolinfillimit):
        infil = pinfil
    else:
        infil = 0.

    return (pinfil - infil)

def surfacepool(wpoolmax, wapl, waps, runon, poolvolremain, poolinfil, delt, fdepth, poolinfillimit, frate, tsurf, lambdaice, rhowater, latentheat, eta):
    """
    Surface Pool

    @param wpoolmax Maximum pool water (liquid plus ice) (mm)
    @param wapl Pool water amount: liquid (mm)
    @param waps Pool water amount: solid (=ice) (mm)
    @param runon Water in excess of what can infiltrate the soil (mm/day)
    @param poolvolremain Unused capacity of surface pool (mm/day)
    @param poolinfil Water flow to pool from other sources than ice thawing (mm/day)    
    @param delt Model time step (day)
    @param fdepth Soil frost layer depth (m)
    @param poolinfillimit Soil frost depth limit for water infiltration (m)
    @param frate Rate of increase of frost layer depth (m/day)
    @param tsurf Soil surface temperature (C)
    @param lambdaice Thermal conductivity of ice (J/m/K/day)
    @param rhowater Density of water (kg/m3)
    @param latentheat Latent heat of water fusion (J/kg)
    @param eta Compound parameter (= LAMBDAice/(RHOwater*LatentHeat)) (m2/K/day)
    @result freezepl Freezing of pool water (mm/day)
    @result thawps Rate of surface ice thawing (mm/day)
    @result pooldrain Water flow from pool to soil (mm/day)
    """
    poolvolremain = max(0., wpoolmax - wapl - waps)
    poolinfil     = min(runon, poolvolremain)
    poolrunoff    = runon - poolinfil
    poolwavail    = poolinfil + wapl/delt
    if (poolwavail == 0.):
        pooldrain = 0.
    elif (fdepth <= poolinfillimit):
        pooldrain = poolwavail
    else :
        pooldrain = max(0.,min( -frate*1000., poolwavail ))

    if ((tsurf>0.) and (wapl==0) and (waps==0.)):
        pirate    = 0.
    else :
        eta       = lambdaice / ( rhowater * latentheat )                                        #! [m2 c-1 day-1]
        pirate    = (sqrt( max(0.,(0.001*waps)**2 - 2.*eta*tsurf*delt)))/delt - (0.001*waps)/delt #! [m day-1]

    if (pirate < 0.):
        freezepl  = 0.
        thawps    = min( waps/delt , -pirate*1000. )
    else :
        freezepl  = max( 0.,min( poolinfil + wapl/delt - pooldrain*delt, pirate*1000. ))
        thawps    = 0.
    return(freezepl, thawps, pooldrain)

def ddayl(doy, lat):
    """
    #!=============================================================================
    #! calculate day length (d d-1) from julian day and latitude (lat, degn)
    #! author - marcel van oijen (ceh-edinburgh)
    #!=============================================================================
    Compute day length

    @param doy Day Of Year (1-366)
    @param lat Latitude (degrees)
    @result dayl Day length (d/d)
    """
    rad  = pi / 180.                                                    #! (radians deg-1)
    dec  = -asin (sin (23.45*rad)*cos (2.*pi*(doy+10.)/365.))           #! (radians)
    decc = max(atan(-1./tan(rad*lat)),min( atan( 1./tan(rad*lat)),dec)) #! (radians)
    dayl = 0.5 * ( 1. + 2. * asin(tan(rad*lat)*tan(decc)) / pi )        #! (d d-1)
    return(dayl)

def penman(dtr, davtmp, wn, vp, lai, rnintc):
    """
    #!=============================================================================
    #! calculate potential rates of evaporation and transpiration (mm d-1)
    #! inputs: lai (m2 m-2), dtr (mj gr m-2 d-1), rnintc (mm d-1)
    #! inputs not in header: vp (kpa), wn (m s-1)
    #! outputs: pevap & ptran (mm d-1)
    #! author - marcel van oijen (ceh-edinburgh)
    #!=============================================================================
    Compute Penman and split E and T

    @param dtr Daily global radiation (MJ GR/m2/day)
    @param davtmp Daily average temperature (Celsius)
    @param wn Windspeed (m/s)
    @param vp Vapour Pressure (kPa)
    @param lai Leaf Area Index (m2 leaf /m2)
    @param rnintc Interception of precipitation by the canopy (mm/day)
    @result pevap Potential rate of evaporation from the soil (mm/day)
    @result ptran Potential transpiration rate (mm/day)
    """
    dtrjm2 = dtr * 1.e6                                    #! (j gr m-2 d-1)
    boltzm = 5.668e-8                                      #! (j m-2 s-1 k-4)
    lhvap  = 2.4e6                                         #! (j kg-1)
    psych  = 0.067                                         #! (kpa degc-1))
    bbrad  = boltzm * pow((davtmp+273.),4) * 86400.        #! (j m-2 d-1)
    svp    = 0.611 * exp(17.4 * davtmp / (davtmp + 239.))  #! (kpa)
    slope  = 4158.6 * svp / pow((davtmp + 239.),2)         #! (kpa degc-1)
    rlwn   = bbrad * max(0.,0.55*(1.-vp/svp))              #! (j m-2 d-1)
    nrads  = dtrjm2 * (1.-0.15) - rlwn                     #! (j m-2 d-1)
    nradc  = dtrjm2 * (1.-0.25) - rlwn                     #! (j m-2 d-1)
    penmrs = nrads * slope/(slope+psych)                   #! (j m-2 d-1)
    penmrc = nradc * slope/(slope+psych)                   #! (j m-2 d-1)
    wdf    = 2.63 * (1.0 + 0.54 * wn)                      #! (kg m-2 d-1 kpa-1)
    penmd  = lhvap * wdf * (svp-vp) * psych/(slope+psych)  #! (j m-2 d-1)
    pevap  =     exp(-0.5*lai)  * (penmrs + penmd) / lhvap #! (mm d-1)
    ptran  = (1.-exp(-0.5*lai)) * (penmrc + penmd) / lhvap #! (mm d-1)
    ptran  = max( 0., ptran-0.5*rnintc )                   #! (mm d-1)
    return(pevap, ptran)

####################################################################################################
# module plant
####################################################################################################

def harvest(rsgcut, tilv, tilg2, claiv, lai, delt, clv, phen, hagere, cres, cst):
    """
    Harvest

    @param rsgcut RS-based observation of grass cut (harvest)
    @param tilv Non-elongating tiller density (-/m2)
    @param tilg2 Elongating generative tiller density (-/m2)
    @param claiv Maxmimum LAI remaining after harvest, when no tillers elongate (m2 leaf / m2)
    @param lai Leaf area index (m2 leaf / m2)
    @param delt Model time step (day)
    @param clv Carbon in leaves (g C/m2)
    @param phen Phenological stage (-)
    @param hagere Fraction of reserves in elongating tillers that is harvested (-)
    @param cres Carbon in reserves (g C/m2)
    @param cst Carbon in stems (g C/m2)
    @result harv Flag indicating that the current day is a harvest day (-)
    @result noharv Flag indicating that the current day is not a harvest day (-)
    @result fractv Fraction of tillers that is not elongating (-)
    @result harvla Harvested leaf area (m2 leaf / m2/day)
    @result harvlv Harvested leaf mass (g C/m2/day)
    @result harvph Resetting of phenological stage by harvesting (-/day)
    @result harvre Harvested reserves (g C/m2/day)
    @result harvst Harvested stem mass (g C/m2/day)
    @result gstub Growth of stubble due to harvest of elongating tillers (g C/m2/day)
    @result harvtilg2 Number of elongating generative tillers harvested today (tiller/day)
    """
    harv   = 0
    noharv = 1
    # Check if we have a harvest today
    if(rsgcut != 0):    
        # We have a harvest !
        harv   = 1
        noharv = 0	
    
    fractv = tilv/(tilg2 + tilv)
    clai   = fractv * claiv
    if (lai <= clai):
        harvfr = 0.0
    else:
        harvfr = 1.0 - clai/lai

    harvla    = (harv   * lai * harvfr) / delt
    harvlv    = (harv   * clv * harvfr) / delt
    harvph    = (harv   * phen        ) / delt
    tv1       = (harvfr * fractv) + (1 - fractv) * hagere
    harvre    = (harv   * tv1 * cres  ) / delt
    harvst    = (harv   * cst         ) / delt
    gstub     = harvst * (1 - hagere)
    harvtilg2 = (harv   * tilg2       ) / delt
    return(harv, noharv, fractv, harvla, harvlv, harvph, harvre, harvst, gstub, harvtilg2)

def biomass(cocresmx, clv, cres, cst):
    """
    Biomass

    @param cocresmx Maximum concentration of reserves in aboveground biomass (not stubble) (-)
    @param clv Carbon in leaves (g C / m2)
    @param cres Carbon in reserves (g C / m2)
    @param cst Carbon in stems (g C / m2)
    @result cresmx Maximum amount of reserves (g C / m2)
    @result resnor Normalised concentration of reserves (-)
    """
    cresmx = cocresmx*(clv + cres + cst)
    resnor = max(0.,min(1., cres/cresmx ))
    return(cresmx, resnor)

def phenology(davtmp, daylp, dayl, daylb, phen, delt, phencr, dlmxge):
    """
    Phenology

    @param davtmp Daily average temperature (Celsius)
    @param daylp Day length below which phenological development slows down (d/d)
    @param dayl Day length (d/d)
    @param daylb Day length below which phenological stage is reset to zero (d/d)
    @param phen Phenological stage (-)
    @param delt Model time step (day)
    @param phencr Phenological stage above which elongation and appearance of leaves on elongating tillers decreases (-)
    @param dlmxge Day length below which DAYLGE becomes less than 1 (d/d)
    @result gphen Rate of phenological development (-/day)
    @result dphen Rate of decrease of phenological stage (-/day)
    @result phenrf Effect of phenological stage on leaf elongation and appearance in elongating tillers (-)
    @result daylge Day length effect on allocation, tillering, leaf appearance, leaf elongation (-)
    """
    gphen = max(0., (davtmp-0.01)*0.000144*24. * (min(daylp,dayl)-0.24) )
    dphen = 0.
    if (dayl < daylb):
        dphen = phen / delt
    phenrf = (1 - phen)/(1 - phencr)
    if (phenrf > 1.0):
        phenrf = 1.0
    if (phenrf < 0.0):
        phenrf = 0.0
    daylge = max(0.,min(1., (dayl - daylb)/(dlmxge - daylb) ))
    return(gphen, dphen, phenrf, daylge)

def foliage1(tbase, davtmp, daylge, slamax, fslamin, resnor):
    """
    Foliage 1

    @param tbase Minimum value of effective temperature for leaf elongation (C)
    @param davtmp Daily average temperature (Celsius)
    @param daylge Day length effect on allocation, tillering, leaf appearance, leaf elongation (-)
    @param slamax Maximum SLA of new leaves (m2 leaf g / C)
    @param fslamin Minimum SLA of new leaves as a fraction of maximum possible SLA (-)
    @param resnor Normalised concentration of reserves (-)
    @result efftmp Effective temperature for leaf elongation (C)
    @result lerv Elongation rate of leaves on non-elongating tillers (m/day)
    @result lerg Elongation rate of leaves on elongating tillers (m/day)
    @result slamin Minimum SLA of new leaves (= SLAMAX * FSLAMIN) (m2 leaf / g C)
    @result slanew SLA of new leaf (m2 leaf / g C)
    """
    efftmp = max(tbase, davtmp)
    lerv   = max(0., (-0.76 + 0.52 * efftmp) / 1000. )
    lerg   = daylge * max(0., (-5.46 + 2.80 * efftmp)/1000. )
    slamin = slamax * fslamin
    slanew = slamax - resnor * (slamax - slamin)
    return(efftmp, lerv, lerg, slamin, slanew)

def lueco2tm(davtmp, rubisc, co2a, kluetilg, fractv, k, parav):
    """
    Calculate luemxq (mol co2 mol-1 par quanta)

    @param davtmp Daily average temperature (Celsius)
    @param rubisc Rubisco content of upper leaves (g/m2 leaf)
    @param co2a CO2 concentration in atmosphere (ppm)
    @param kluetilg LUE-increase with increasing fraction elongating tillers (-)
    @param fractv Fraction of tillers that is not elongating (-)
    @param k PAR extinction coefficient (m2 / m2 leaf)
    @param parav Average PAR during the photoperiod (mumol PAR /m2/s) (micromol par quanta m-2 s-)
    @result rehardredstart
    @result resphardsi
    """
    t      = davtmp                                            #!(degc)
    rubiscn = rubisc * (1.e6/550000.)
    eavcmx =  68000                                            #!(j mol-1)
    eakmc  =  65800                                            #!(j mol-1)
    eakmo  =   1400                                            #!(j mol-1)
    kc25   =     20                                            #!(mol co2 mol-1 rubisco s-1)
    kmc25  =    460                                            #!(ppm co2)
    kmo25  =     33                                            #!(% o2)
    kokc   =      0.21                                         #!(-)
    o2     =     21                                            #!(% o2)
    r      =      8.314                                        #!(j k-1 mol-1)
    co2i   = 0.7 * co2a                                        #!(ppm co2)
    vcmax  = rubiscn * kc25 * exp((1/298.-1/(t+273))*eavcmx/r) #!(micromol co2 m-2 s-1)
    kmc    =         kmc25 * exp((1/298.-1/(t+273))*eakmc /r)  #!(ppm co2)
    kmo    =         kmo25 * exp((1/298.-1/(t+273))*eakmo /r)  #!(% o2)
    gammax = 0.5 * kokc * kmc * o2 / kmo                       #!(ppm co2)
    pmax   = vcmax * (co2i-gammax) / (co2i + kmc * (1+o2/kmo)) #!(micromol co2 m-2 s-1)
    tmpfac = max( 0., min( 1., (t+4.)/5. ) )                   #!(-)
    eff    = tmpfac * (1/2.1) * (co2i-gammax) / (4.5*co2i+10.5*gammax) #!(mol co2 mol-1 par quanta)
    luemxq = eff*pmax*(1+kluetilg*(1-fractv)) / (eff*k*parav + pmax) #!(mol co2 mol-1 par quanta)aug 8    
    return(luemxq)
  
def hardeningsink(rehardredend, rehardredday, doy, tsurf, thardmx, lt50, lt50mn, hparam, clv, kresphard, resnor):
    rehardredstart = (rehardredend-rehardredday)%365. #Modulo operation
    doysincestart  = (doy-rehardredstart)%365. #Modulo operation
    if ( doysincestart < (rehardredday + 0.5 * (365. - rehardredday)) ):
        rehardperiod = max( 0., 1. - doysincestart/rehardredday)
    else:
        rehardperiod = 1.

    if ( (tsurf > thardmx) or (lt50 < lt50mn) ):
        rateh = 0.
    else:
        rateh = rehardperiod * hparam * (thardmx - tsurf) * (lt50 - lt50mn)

    resphardsi = rateh * clv * kresphard * max(0.,min(1., resnor*5.))
    return(rehardredstart, resphardsi)


def growth(parint, tranrf, luemxq, noharv, cres, tcres, davtmp, resphardsi, cresmx, delt, tilg2, cst, cstavm, simax1t, phenrf, nellvm, lerv, tilv, lfwidv, lerg, lfwidg, shape, slanew, yg, daylge):
    """
    Growth

    @param parint PAR interception (mol PAR /m2/day)
    @param tranrf Transpiration realisation factor (-)
    @param luemxq Light-use efficiency (mol CO2/mol PAR)
    @param noharv Flag indicating that the current day is not a harvest day (-)
    @param cres Carbon in reserves (g C/m2)
    @param tcres Time constant of mobilisation of reserves (day)
    @param davtmp Daily average temperature (Celsius)
    @param resphardsi Sink strength from carbohydrate demand of hardening (g C/m2/day)
    @param cresmx Maximum amount of reserves (g C/m2)
    @param delt Model time step (day)
    @param tilg2 Elongating generative tiller density (-/m2)
    @param cst Weight of stems (g C/m2)
    @param cstavm Maximum size of elongating tillers (g C/tiller)
    @param simax1t Sink strength of small elongating tillers (g C/tiller)
    @param phenrf Effect of phenological stage on leaf elongation and appearance in elongating tillers (-)
    @parma nellvm Number of elongating leaves per non-elongating tiller (-/tiller)
    @param lerv Elongation rate of leaves on non-elongating tillers (m/day)
    @param tilv Non-elongating tiller density (-/m2)
    @param lfwidv Leaf width on non-elongating tillers (m)
    @param lerg Elongation rate of leaves on elongating tillers (m/day)
    @parma lfwidg Leaf width on elongating tillers (m)
    @param shape Area of a leaf relative to a rectangle of same length and width (-)
    @param slanew SLA of new leaf (m2 leaf g/C-)
    @param yg Growth yield per unit expended carbohydrate (g C / g C)
    @param daylge Day length effect on allocation, tillering, leaf appearance, leaf elongation (-)
    @result resmob Mobilisation of reserves (g C/m2/day)
    @result resphard Plant hardening respiration (g C/m2/day)
    @result gres Gross growth rate of reserve pool, uncorrected for remobilisation (g C/m2/day)
    @result glv Growth of leaf mass (g C/m2/day)
    @result gst Growth of stems (g C/m2/day)
    @result grt Growth of roots (g C/m2/day)
    @result respgsh Respiration associated with shoot growth (g C/m2/day)
    @result respgrt Respiration associated with root growth (g C/m2/day)
    """
    phot     = parint * tranrf * 12. * luemxq * noharv
    resmob   = (cres * noharv / tcres) * max(0., min( 1., davtmp/5. ))
    source   = resmob + phot
    resphard = min(source, resphardsi)
    allotot  = source - resphard
    gressi   = 0.5 * (resmob + max(0., (cresmx - cres) / delt))
    if (tilg2 != 0.0):
        cstav  = cst / tilg2 
    else:
        cstav  = 0.
    sink1t   = max(0., 1 - (cstav / cstavm)) * simax1t
    nellvg   = phenrf * nellvm 
    glaisi   = ((lerv * tilv * nellvm * lfwidv) + (lerg * tilg2 * nellvg * lfwidg)) * shape * tranrf
    glvsi    = (glaisi * noharv / slanew) / yg
    gstsi    = (sink1t * tilg2 * tranrf * noharv) / yg
    # Run allocation
    gres, glv, gst, grt, respgsh, respgrt = allocation(glvsi, gstsi, daylge, allotot, gressi, yg)
    return(nellvg, resmob, resphard, gres, glv, gst, grt, respgsh, respgrt)


def allocation(glvsi, gstsi, daylge, allotot, gressi, yg):
    """
    Allocation

    @param glvsi Potential growth rate of leaf mass (g C/m2/day)
    @param gstsi Potential growth rate of stems (g C/m2/day)
    @param daylge Day length effect on allocation, tillering, leaf appearance, leaf elongation (-)
    @param allotot Allocation of carbohydrates to sinks other than hardening (g C/m2/day)
    @param gressi Sink strength of reserve pool (g C/m2/day)
    @param yg Growth yield per unit expended carbohydrate (g C/G C)
    @result gres Gross growth rate of reserve pool, uncorrected for remobilisation (g C/m2/day)
    @result glv Growth of leaf mass (g C/m2/day)
    @result gst Growth of stems (g C/m2/day)
    @result grt Growth of roots (g C/m2/day)
    @result respgsh Respiration associated with shoot growth (g C/m2/day)
    @result respgrt Respiration associated with root growth (g C/m2/day)
    """
    gshsi = glvsi + gstsi
    if (daylge >= 0.1):
        #! situation 1: growth has priority over storage (spring and growth period)
        #! calculate amount of assimilates allocated to shoot
        allosh = min( allotot, gshsi )
        #! calculate amount of assimilates allocated to reserves    
        gres   = min( allotot - allosh, gressi)
    else:
        #! situation 2: storage has priority over shoot (autumn)
        #! calculate amount of assimilates allocated to reserves
        gres   = min( allotot, gressi )
        #! calculate amount of assimilates allocated to shoot
        allosh = min( allotot - gres, gshsi )
     
    #! all surplus carbohydrate goes to roots
    allort  = allotot - allosh - gres
    if (gshsi == 0.):
        gshsi = 1

    allolv  = glvsi * (allosh / gshsi)
    allost  = gstsi * (allosh / gshsi)
    glv     = allolv * yg
    gst     = allost * yg
    grt     = allort * yg
    respgsh = (allolv + allost) * (1-yg)
    respgrt =  allort           * (1-yg)
    return(gres, glv, gst, grt, respgsh, respgrt)
    
def plantrespiration(fo2, fo2mx, respgrt, respgsh, resphard):
    """
    Plant Respiration

    @param fo2 Soil oxygen as a fraction of total gas (mol O2/mol gas)
    @param fo2mx Maximum oxygen fraction of soil gas (mol O2/mol gas)
    @param respgrt Respiration associated with root growth (g C/m2/day)
    @param respgsh Respiration associated with shoot growth (g C/m2/day)
    @param resphard Plant hardening respiration (g C/m2/day)    
    @result rplantaer Aerobic plant respiration (g C/m2/day)
    """
    faer      = max(0.,min(1., fo2/fo2mx ))
    rplantaer = faer * ( respgrt + respgsh + resphard )
    return(rplantaer)

def senescence(permgas, tanaer, delt, ldt50a, ldt50b, lt50, krdranaer, krsr3h, tsurf, dparam, lt50mx, tsurfdiff, ratedmx, resphard, clv, kresphard, lai,laicr, rdrsco, rdrsmx, rdrtmin, rdrtem, noharv, cstub, rdrstub, tilv, crt, rdroot):
    """
    Senescence Function

    @param permgas Permeability of soil surface to gas exchange (-/day)
    @param tanaer
    @param delt Model time step (day)
    @param ldt50a Intercept of linear dependence of LD50 on lT50 (day)
    @param ldt50b Slope of linear dependence of LD50 on LT50 (day/C)
    @param lt50 Sensitivity to Frost (Lethal Temperature 50%) Temperature that kills half the plants in a day (C)
    @param krdranaer Maximum relative death rate due to anearobic conditions (-/day)
    @param krsr3h Constant in the logistic curve for frost survival (-/C)
    @param tsurf Soil surface temperature (C)
    @param dparam Constant in the calculation of dehardening rate (-/C/day)
    @param lt50mx Maximum LT50 (C)
    @param tsurf Soil surface temperature (C)diff
    @param ratedmx Maximum dehardening rate (C/day)
    @param resphard Plant hardening respiration (g C/m2/day)    
    @param clv Carbon in leaves (g C/m2)
    @param kresphard Carbohydrate requirement of hardening (g C/ g C /C)
    @param lai Leaf Area Index (m2 leaf / m2)
    @param laicr LAI above which shading induces leaf senescence (m2 leaf / m2)
    @param rdrsco Relative death rate of leaves and non-elongating tillers due to shading when LAI is twice the threshold (LAICR) (-/day)
    @param rdrsmx Maximum relative death rate of leaves and non-elongating tillers due to shading (-/day)
    @param rdrtmin Minimum relative death rate of foliage (-/day)
    @param rdrtem Proportionality of leaf senescence with temperature (-/day/C)
    @param noharv Flag indicating that the current day is not a harvest day (-)
    @param cstub Carbon in stubble (g C/m2)
    @param rdrstub Relative death rate of stubble (-/day)
    @param tilv Non-elongating tiller density (-/m2)
    @param crt Carbon in roots (g C/m2)
    @param rdroot Relative death rate of roots (-/day)
    @result dtanaer Change in days since start anaerobic conditions (day/day)
    @result hardrate Hardening (LT50 becoming more negative) (C/day)
    @result dehardrate Dehardening rate (LT50 becoming less negative) (C/day)
    @result tv2til Relative death rate of non-elongating tillers (-/day)
    @result dlai Death rate of leaf area (m2 leaf /m2/day)
    @result dlv Death rate of leaf mass (g C/m2/day)
    @result dstub Death rate of stubble (g C/m2/day)
    @result dtilv Death rate of non-elongating tillers (tillers/m2/day)
    @result drt Death rate of roots (g C/m2/day)
    """
    dtanaer, rdrtox = anaerobicdamage(permgas, tanaer, delt, ldt50a, ldt50b, lt50, krdranaer)
    rdrfrost, hardrate, dehardrate = hardening(krsr3h, tsurf, lt50, dparam, lt50mx, tsurfdiff, delt, ratedmx, resphard, clv, kresphard)
    if (lai < laicr):
        tv1 = 0.0 
    else:
        tv1 = rdrsco * (lai - laicr) / laicr

    rdrs   = min(tv1, rdrsmx)
    rdrt   = max(rdrtmin, rdrtem * tsurf)
    tv2    = noharv * max(rdrs, rdrt, rdrfrost, rdrtox)
    tv2til = noharv * max(rdrs, rdrfrost, rdrtox)
    dlai   = lai    * tv2
    dlv    = clv    * tv2
    dstub  = cstub  * rdrstub
    dtilv  = tilv   * tv2til
    drt    = crt    * rdroot
    return(dtanaer, hardrate, dehardrate, tv2til, dlai, dlv, dstub, dtilv, drt)

def anaerobicdamage(permgas, tanaer, delt, ldt50a, ldt50b, lt50, krdranaer):
    """
    Anaerobic Damage

    @param permgas Permeability of soil surface to gas exchange (-/day)
    @param tanaer Time since start anaerobic conditions (d)
    @param delt Model time step (day)
    @param ldt50a Intercept of linear dependence of LD50 on lT50 (day)
    @param ldt50b Slope of linear dependence of LD50 on LT50 (day/C)
    @param lt50 Sensitivity to Frost (Lethal Temperature 50%) Temperature that kills half the plants in a day (C)
    @parma krdranaer Maximum relative death rate due to anearobic conditions (-/day)
    @result dtanaer Change in days since start anaerobic conditions (day/day)
    @result rdrtox Relative death rate due to anaerobic conditions (day/day)
    """
    if (permgas == 0.):
        dtanaer = 1.
    else:
        dtanaer = -tanaer / delt
    
    ld50 = ldt50a + ldt50b * lt50
    if (tanaer > 0.):
        rdrtox = krdranaer / (1. + exp( -krdranaer * (tanaer - ld50)))
    else:
        rdrtox = 0.

    return(dtanaer, rdrtox)
     
def hardening(krsr3h, tsurf, lt50, dparam, lt50mx, tsurfdiff, delt, ratedmx, resphard, clv, kresphard):
    """
    Hardening

    @param krsr3h Constant in the logistic curve for frost survival (-/C)
    @param tsurf Soil surface temperature (C)
    @param lt50 Sensitivity to Frost (Lethal Temperature 50%) Temperature that kills half the plants in a day (C)
    @param dparam Constant in the calculation of dehardening rate (-/C/day)
    @param lt50mx Maximum LT50 (C)
    @param tsurf Soil surface temperature (C)diff
    @param delt Model time step (day)
    @param ratedmx Maximum dehardening rate (C/day)
    @param resphard Plant hardening respiration (g C/m2/day)    
    @param clv Carbon in leaves (g C/m2)
    @param kresphard Carbohydrate requirement of hardening (g C/ g C /C)
    @result rdrfrost Relative death rate due to frost (-/day)
    @result hardrate Hardening (LT50 becoming more negative) (C/day)
    @result dehardrate Dehardening rate (LT50 becoming less negative) (C/day)
    """
    rsr3h      = 1. / (1. + exp( -krsr3h * (tsurf-lt50) ))
    #! rdrfrost should be less than 1 to avoid numerical problems
    #! (loss of all biomass but keeping positive reserves). we cap it at 0.5.
    rsrday     = rsr3h #! in previous versions we had rsrday = rsr3h^8 which underestimated survival
    rdrfrost   = min( 0.5, 1. - rsrday )
    rated      = min( dparam * (lt50mx - lt50) * (tsurf + tsurfdiff), (lt50mx - lt50) / delt )
    dehardrate = max(0.,min( ratedmx, rated ))
    hardrate   = resphard / (clv * kresphard)
    return(rdrfrost, hardrate, dehardrate)

def foliage2(slanew, glv, tsurf, tbase, phy, noharv, tranrf, daylge, fractv, phenrf, laitil, laieft, lai, fsmax, resnor, tilv, davtmp, toptge, tv2til, rgenmx, vern, dayl, daylg1g2, tilg1, rgrtg1g2):
    """
    Foliage 2

    @param slanew SLA of new leaf (m2 leaf g/C-)
    @param glv Growth rate of leaf mass (g C/m2/day)
    @param tsurf Soil surface temperature (C)
    @param tbase Minimum value of effective temperature for leaf elongation (C)
    @param phy Phyllochron (C day)
    @param noharv Flag indicating that the current day is not a harvest day (-)
    @param tranrf Transpiration realisation factor (-)
    @param daylge Day length effect on allocation, tillering, leaf appearance, leaf elongation (-)
    @param fractv Fraction of tillers that is not elongating (-)
    @param phenrf Effect of phenological stage on leaf elongation and appearance in elongating tillers (-)
    @param laitil Maximum ratio of tiller and leaf apearance at low leaf area index (-)
    @param laieft Decrease in tillering with leaf area index (m2/m2 leaf)
    @param lai Leaf Area Index (m2 leaf / m2)
    @param fsmax Maximum ratio of tiller and leaf appearance based on sward geometry (-)
    @param resnor Normalised concentration of reserves (-)
    @param tilv Non-elongating tiller density (-/m2)
    @param davtmp Daily average temperature (Celsius)
    @param toptge Optimum temperature for vegetative tillers becoming generative (C)
    @param tv2til Relative death rate of non-elongating tillers (-/day)
    @param rgenmx Maximum relative rate of tillers becoming elongating tillers (-/day)
    @param vern Vernalization status (-)
    @param dayl Day length (day/day)
    @param daylg1g2 Minimum day length above which generative tillers can start elongating (by moving from TILG1 to TILG2). (day/day)
    @param tilg1 Non-elongating generative tiller density (-/m2)
    @param rgrtg1g2 Relative rate of TILG1 becoming TILG2 (tiller/tiller/day)
    @result glai Growth rate of leaf area (m2 leaf /m2/day)
    @result gtilv Growth of non-elongating tillers (g C/m2/day)
    @result tilvg1 Number of non-elongating tillers (tillers)
    @result tilg1g2 Number of Non-elongating generative tiller becoming elongating generative tillers (tillers)
    """
    glai    = slanew * glv
    if (tsurf < tbase):
        tv1   = 0.
    else:
        tv1   = tsurf/phy
        
    rleaf   = tv1 * noharv * tranrf * daylge * ( fractv + phenrf * (1 - fractv) )
    tv2     = laitil - laieft * lai
    if (tv2 > fsmax):
        tv2 = fsmax

    if (tv2 < 0.):
        tv2 = 0.

    rgrtv   = max( 0, tv2 * resnor * rleaf )
    gtilv   = tilv  * rgrtv
    tge     = max( 0, 1 - (abs(davtmp - toptge)) / (toptge - tbase))
    rgrtvg1 = min( 1 - tv2til, noharv * daylge * tge * rgenmx ) * vern
    tilvg1  = tilv  * rgrtvg1
    if (dayl > daylg1g2):
        tilg1g2 = tilg1 * rgrtg1g2
    else:
        tilg1g2 = 0.
    return(glai, gtilv, tilvg1, tilg1g2)

####################################################################################################
# Soil Module
####################################################################################################

def soilwatercontent(fdepth, rootd, wal):
    """
    Estimate and return wcl

    @param fdepth Soil frost layer depth (m)
    @param rootd Rooting depth (m)
    @param wal Soil water amount: liquid (mm)
    @result wcl Water concentration in non-frozen soil (m3/m3)
    """
    if (fdepth < rootd):
        return wal * 0.001 / (rootd-fdepth)
    else:
        return 0

def physics(davtmp, fdepth, rootd, sdepth, was, gamma, wcl, delt, wcfc, tsurf, lambdasoil, rhowater, latentheat):
    """
    Physics update

    @param davtmp Daily average temperature (Celsius)
    @param fdepth Soil frost layer depth (m)
    @param rootd Rooting depth (m)
    @param sdepth Snow depth (m)
    @param was Soil water amount: solid (=ice) (mm)
    @param gamma Temperature extinction coefficient of snow (-/m)
    @param wcl Water concentration in non-frozen soil (m3/m3)
    @param delt Model time step (day)
    @param wcfc Water concentration at field capacity (m3/m3)
    @param tsurf Soil surface temperature (C)
    @param lambdasoil Soil parameter (J/m/C/day)
    @param rhowater Density of water (kg/m3)
    @param latentheat Latent heat of water fusion (J/kg)
    @result tsurf Soil surface temperature (C)
    @result fperm Flag for full soil depth Freezing (0/1)
    @result frate Rate of increase of frost layer depth (m/day)
    """
    # update values of tsurf, fperm, and frate (frozen soil function call)
    if (fdepth > 0.):
        tsurf = davtmp / (1. + 10. * (sdepth / fdepth) )
        fperm = 0.
    else :
        tsurf = davtmp * exp(-gamma*sdepth)
        fperm = 1.
   
    frate = frozensoil(fdepth, rootd, wcfc, was, wcl, tsurf, lambdasoil, rhowater, latentheat, delt)
    return(tsurf, fperm, frate)

  
def frozensoil(fdepth, rootd, wcfc, was, wcl, tsurf, lambdasoil, rhowater, latentheat, delt):
    """
    Frozen Soil Variation

    @param fdepth Soil frost layer depth (m)
    @param rootd Rooting depth (m)
    @param wcfc Water concentration at field capacity (m3/m3)
    @param was Soil water amount: solid (=ice) (mm)
    @param wcl Water concentration in non-frozen soil (m3/m3)
    @param tsurf Soil surface temperature (C)
    @param lambdasoil Soil parameter (J/m/C/day)
    @param rhowater Density of water (kg/m3)
    @param latentheat Latent heat of water fusion (J/kg)
    @param delt Model time step (day)
    @result Frozem Soil variation
    """
    #! determining the amount of solid water that contributes in transportation of heat to surface 'wceff'
    if (fdepth > rootd):
         wceff = wcfc
    elif (fdepth > 0.):
         wceff = (0.001*was) / fdepth
    else:
         wceff = wcl
    #! calculating potential frost rate 'pfrate'
    if (((fdepth == 0.) and (tsurf>0.)) or (wceff == 0.)): #! no soil frost present and no frost starting
       pfrate = 0.
    else:
       alpha  = lambdasoil / ( rhowater * wceff * latentheat )
       pfrate = sqrt( max(0.,fdepth**2 - 2.*alpha*tsurf) ) - fdepth
     
    if ((pfrate >= 0.) and (fdepth > 0.) and (fdepth < rootd)):
       return(pfrate * (0.001*was/fdepth) / wcfc) # soil frost increasing
    elif ((pfrate+fdepth/delt) < 0.):
       return(-fdepth / delt)                     # remaining soil frost thaws away
    else :
       return(pfrate)
     
def frdrunir(wcfc, rootd, fdepth, wcst, infil, pooldrain, wal, delt, evap, tran, frate, was, drate, irrigf):
    """
    Freezing, Thawing, Draining, Runoff, Irrigation

    @param wcfc Water concentration at field capacity (m3/m3)
    @param rootd Rooting depth (m) 
    @param fdepth Soil frost layer depth (m)
    @param wcst Water concentration at saturation (m3/m3)
    @param infil Water flow into soil from precipitation and snow melt (mm/day)
    @param pooldrain Water flow from pool to soil (mm/day)
    @param wal Soil water amount: liquid (mm)
    @param delt Model time step (day)
    @param evap Evaporation of water from soil surface (mm/day)
    @param tran Transpiration of water from soil surface (mm/day)
    @param frate Rate of increase of frost layer depth (m/day)
    @param was Soil water amount: solid (=ice) (mm)
    @param drate Maximum drainage rate (mm/day)
    @param irrigf Irrigation relative to what would maintain field capacity (-)
    @result freezel Freezing of soil water (mm/day)
    @result thaws Water flow to soil pool from thawing of frozen soil (mm/day)
    @result drain Drainage rate below the root zone (mm/day)
    @result runoff Loss of water by runoff (mm/day)
    @result irrig Irrigation rate (mm/day)
    """
    wafc   = 1000. * wcfc * max(0.,(rootd-fdepth))                     # ! (mm)
    wast   = 1000. * wcst * max(0.,(rootd-fdepth))                     # ! (mm)
    infiltot = infil + pooldrain
    if (fdepth < rootd):
        freezel = max(0., min( wal/delt + (infiltot - evap - tran), (frate/(rootd-fdepth))*wal)) # ! (mm d-1)
    else:
        freezel = 0.                                                     # ! (mm d-1)
  
    if ((fdepth > 0.) and (fdepth <= rootd)):
        thaws   = max(0.,min( was/delt, -frate*was/fdepth ))             # ! (mm d-1)
    else :
        thaws   = 0.                                                     # ! (mm d-1)
  
    drain  = max(0.,min( drate, (wal-wafc)/delt + (infiltot - evap - tran - freezel + thaws) )) #! (mm d-1)
    runoff = max(0.,(wal-wast)/delt + (infiltot - evap - tran - freezel + thaws - drain) ) #! (mm d-1)
    irrig  = irrigf * ((wafc-wal)/delt - (infiltot - evap - tran - freezel + thaws - drain - runoff))#! (mm d-1)
    return(freezel, thaws, drain, runoff, irrig)

def o2status(o2, rootd, fgas):
    """
    O2 Status

    @param o2 Soil oxygen content (mol/m2)
    @param rootd Rooting depth (m)
    @param fgas Soil pore space (potentially gaseous) (m3/m3)
    @result O2 Status
    """
    return( o2 / (rootd * fgas * 1000./22.4))
  
def o2fluxes(fo2mx, rootd, fgas, rplantaer, krtotaer, permgas, o2, delt):
    """
    O2 Fluxes

    @param fo2mx Maximum oxygen fraction of soil gas (mol O2/mol gas)
    @param rootd Rooting depth (m)
    @param fgas Soil pore space (potentially gaseous) (m3/m3)
    @param rplantaer Aerobic plant respiration (g C/m2/day)
    @param krtotaer Ratio of total to aerobic respiration (-)
    @param permgas Permeability of soil surface to gas exchange (-/day)
    @param o2 Soil oxygen content (mol/m2)
    @param delt Model time step (day)
    @result o2mx Maximum oxygen content of soil (mol/m2)
    @result o2out Efflux of oxygen from the soil (mol/m2/day)
    @result o2in Influx of oxygen into the soil (mol/m2/day)
    """
    o2mx  = fo2mx * rootd * fgas * 1000./22.4
    o2out = rplantaer * krtotaer * 1./12. * 1.
    o2in =  permgas * ( (o2mx - o2) + o2out * delt )
    return(o2mx, o2out, o2in)

####################################################################################################
# Resources Module
####################################################################################################

def light(dayl, par, k, lai, dtr):
    """
    Light

    @param dayl Day Length (day/day)
    @param par Daily photosynthetically active radiation (mol PAR/m2/day)
    @param k PAR extinction coefficient (m2 / m2 leaf)
    @param lai Leaf Area index (m2 leaf / m2)
    @param dtr Daily global radiation (MJ GR/m2/day)
    @result parav Average PAR during the photoperiod (mumol PAR /m2/s) (micromol par quanta m-2 s-)
    @result parint PAR interception (mol PAR /m2/day)
    @result dtrint Interception of global radiation (MJ GR/m2/day)
    """
    if (dayl > 0):
        parav = par * (1e6/(24*3600)) / dayl
    else :
        parav = 0.
    
    parint = par * (1 - exp(-1.0*k*lai))
    dtrint = dtr * (1 - exp(-0.75*k*lai))
    return(parav, parint, dtrint)

def evaptrtrf(fdepth, rootd, wal, wcad, pevap, wcfc, wcwp, ptran, tranco, wcst, wcwet, delt):
    """
    Evaporation, Transpiration & Transpiration Fraction

    @param fdepth Soil frost layer depth (m)
    @param rootd Rooting depth (m)
    @param wal Soil water amount: liquid (mm)
    @param wcad Water concentration at air dryness (m3/m3)
    @param pevap Potential rate of evaporation from the soil (mm/day)
    @param wcfc Water concentration at field capacity (m3/m3)
    @param wcwp Water concentration at wilting point (m3/m3)
    @param ptran Potential transpiration rate (mm/day)
    @param tranco Transpiration coefficient (mm/day)
    @param wcst Water concentration at saturation (m3/m3)
    @param wcwet Water concentration above which transpiration is reduced (m3/m3)
    @param delt Model time step (day)
    @result evap Evaporation of water from soil surface (mm/day)
    @result tran Transpiration of water from soil surface (mm/day)
    @result tranrf Transpiration realisation factor (-)
    """
    if (fdepth < rootd):
        wcl = wal*0.001 / (rootd-fdepth)
    else:
        wcl = 0                                                 # ! (m3 m-3)
    waad = 1000. * wcad * (rootd-fdepth)                        # ! (mm)
    evap = pevap * max(0., min(1., (wcl-wcad)/(wcfc-wcad) ))    # ! (mm d-1)
    wccr = wcwp + max( 0.01, ptran/(ptran+tranco)*(wcfc-wcwp))  #! (m3 m-3)
    if (wcl > wccr):
        if(wcst == wcwet):
            fr = 1
        else:
            fr = max(0., min(1., (wcst-wcl)/(wcst-wcwet)))              
    else:
        if(wccr == wcwp):
            fr = 1
        else:
            fr = max(0., min(1., (wcl-wcwp)/(wccr-wcwp)))       #! (mm mm-1)

    tran = ptran * fr                                           #! (mm d-1)
    if (evap+tran > 0.):
        availf = min( 1., ((wal-waad)/delt) / (evap+tran) )         
    else:
        availf = 0                                              #! (mm mm-1)
    
    evap = evap * availf                                        #! (mm d-1)
    tran = tran * availf                                        #! (mm d-1)
    if (ptran > 0.):
        tranrf = tran / ptran                                   #! (-)
    else:
        tranrf = 1                                              #! (-)
    
    return(evap, tran, tranrf)

def rootdg(fdepth, rootd, wal, rootdm, wcwp, rrdmax, delt, wcfc):
    """
    Root Depth "g?" 
    
    @param fdepth Soil frost layer depth (m)
    @param rootd Rooting depth (m)
    @param wal Soil water amount: liquid (mm)
    @param rootd Rooting depth (m)m Initial and maximum value rooting depth (m)
    @param wcwp Water concentration at wilting point (m3/m3)
    @param rrdmax Maximum root depth growth rate (m/day)
    @param delt Model time step (day)
    @param wcfc Water concentration at field capacity (m3/m3)
    @result wcl Water concentration in non-frozen soil (m3/m3)
    @result rrootd Root depth growth rate (m/day)
    @result explor Increased access to water by root depth growth (m/day)
    """
    if (fdepth < rootd):
        wcl = wal*0.001 / (rootd-fdepth)
    else:
        wcl = 0                                                 #! (m3 m-3)
    
    if ((rootd < rootdm) and (wcl > wcwp)):
        rrootd = min( rrdmax, (rootdm-rootd)/delt )
    else:
        rrootd = 0.
    
    explor = 1000. * rrootd * wcfc
    return(wcl, rrootd, explor)


