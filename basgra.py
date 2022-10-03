from lib_basgra import *

def basgra(params, weather, startdoy, enddoy, startyear, endyear):
    """
    BASGRA model as a function
    #-------------------------------------------------------------------------------
    # this is the basic grass model originally written in matlab/simulink by marcel
    # van oijen, mats hoglind, stig morten thorsen and ad schapendonk.
    # 2011-07-13: translation to fortran by david cameron and marcel van oijen.
    # 2014-03-17: extra category of tillers added
    # 2014-04-03: vernalization added
    # 2014-04-03: lower limit of temperature-driven leaf senescence no longer zero
    #-------------------------------------------------------------------------------
    @param params Parameters for BASGRA configration
    @param weather Weather file (with harvest days and optional RS input)
    @result harvestedBiomass The biomass harvested estimation of the model
    """
    # initial constants
    log10clvi   = params[0]
    log10cresi  = params[1]
    log10crti   = params[2]
    csti        = params[3]
    log10laii   = params[4]	   
    pheni       = params[5]
    tiltoti     = params[6]
    frtilgi     = params[7]
    lt50i       = params[8]
    #process parameters 
    claiv       = params[9]	   
    cocresmx    = params[10]
    cstavm      = params[11]
    daylb       = params[12] 
    daylp       = params[13]
    dlmxge      = params[14]
    fslamin     = params[15]	 
    fsmax       = params[16]
    hagere      = params[17]
    k           = params[18]
    laicr       = params[19]
    laieft      = params[20]   
    laitil      = params[21]
    lfwidg      = params[22]
    lfwidv      = params[23]
    nellvm      = params[24]	 
    phencr      = params[25]	  
    phy         = params[26]	  
    rdrsco      = params[27]  
    rdrsmx      = params[28]	 
    rdrtem      = params[29]
    rgenmx      = params[30] 
    rootdm      = params[31]	  
    rrdmax      = params[32]
    rubisc      = params[33]	   
    shape       = params[34]	  
    simax1t     = params[35]
    slamax      = params[36]
    tbase       = params[37]	   
    tcres       = params[38]	  
    toptge      = params[39]	  
    tranco      = params[40]	 
    yg          = params[41]
    lat         = params[42]
    wci         = params[43]
    fwcad       = params[44]
    fwcwp       = params[45]
    fwcfc       = params[46]
    fwcwet      = params[47]
    wcst        = params[48]
    wpoolmax    = params[49]
    dparam      = params[50]		
    fgas        = params[51]		
    fo2mx       = params[52]		
    gamma       = params[53]
    hparam      = params[54]		
    krdranaer   = params[55]		
    kresphard   = params[56]		
    krsr3h      = params[57]		
    krtotaer    = params[58]		
    ksnow       = params[59]		
    lambdasoil  = params[60]		
    ldt50a      = params[61]		
    ldt50b      = params[62]		
    lt50mn      = params[63]		
    lt50mx      = params[64]		
    ratedmx     = params[65]		
    rehardredday= params[66]		
    rhonewsnow  = params[67]			
    rhopack     = params[68]		
    swret       = params[69]		
    swrf        = params[70]		
    thardmx     = params[71]		
    tmeltfreeze = params[72]			
    trainsnow   = params[73]			
    tsurfdiff   = params[74]
    kluetilg    = params[75]
    frtilgg1i   = params[76]
    daylg1g2    = params[77]
    rgrtg1g2    = params[78]
    rdrtmin     = params[79]
    tvern       = params[80]
    clitt0      = params[81] #300
    csom0       = params[82] #10000
    cnlitt0     = params[83] #50
    cnsomf0     = params[84] #25
    cnsoms0     = params[85] #25
    fcsomf0     = params[86] #0.5
    flittsomf   = params[87] #0.6
    fsomfsoms   = params[88] #0.03
    rnleach     = params[89] #0.5
    knemit      = params[90] #0.0005
    nmin0       = params[91] #1
    tclitt      = params[92] #365
    tcsomf      = params[93] #1460
    tcsoms      = params[94] #73000
    tmaxf       = params[95] #50
    tsigmaf     = params[96] #20
    rfn2o       = params[97] #8.5
    wfps50n2o   = params[98] #0.7
    ncshmax     = params[99] #0.04
    ncr         = params[100] #0.03
    rdroot      = params[101] #0.02
    rdrstub     = params[102] #0.2
    nfertmult   = params[103] #1
    fncgshmin   = params[104] #0
    tcnshmob    = params[105] #8
    tcnupt      = params[106] #8
    f_digest_wall_fmin  = params[107] #0.7
    f_digest_wall_max   = params[108] #0.9
    f_wall_lv_fmin      = params[109] #0.7
    f_wall_lv_max       = params[110] #0.6
    f_wall_st_fmin      = params[111] #0.8

    # parameter transformations 
    clvi  = pow(10,log10clvi)
    cresi = pow(10,log10cresi)
    crti  = pow(10,log10crti)
    laii  = pow(10,log10laii)
    wcad  = fwcad  * wcst
    wcwp  = fwcwp  * wcst
    wcfc  = fwcfc  * wcst
    wcwet = fwcwet * wcst
    # End of input parameters

    # calendar & weather
    yeari   = weather[:,0]  # YYYY
    doyi    = weather[:,1]  # DOY
    gri     = weather[:,2]  # Global radiation Incoming
    tmmni   = weather[:,3]  # Tmin
    tmmxi   = weather[:,4]  # Tmax
    vpi     = weather[:,5]  # Vapour Pressure
    raini   = weather[:,6]  # Rainfall
    wni     = weather[:,7]  # Wind Speed
    rsevapi = weather[:,8]  # RS evaporation [ optional, zero values if not used ]
    rstrani = weather[:,9]  # RS transpiration [ optional, zero values if not used ]
    rslaii  = weather[:,10] # RS LAI [ optional, zero values if not used ]
    rsgcuti = weather[:,11] # RS grass cut

    # initial constants
    clv     = clvi
    clvd    = 0             # was clvdi=0 in parameters_plant.f90
    yielD   = 0             # Renamed yielD bc yield is a Python built-in function (reserved)
    yield_last = 0
    cres    = cresi
    crt     = crti
    cst     = csti
    cstub   = 0
    drystor = 0
    fdepth  = 0
    lai     = laii
    lt50    = lt50i
    o2      = fgas * rootdm * fo2mx * 1000./22.4
    phen    = pheni
    rootd   = rootdm
    sdepth  = 0
    tanaer  = 0
    tilg1   = tiltoti * frtilgi * frtilgg1i
    tilg2   = tiltoti * frtilgi * (1-frtilgg1i)
    tilv    = tiltoti * (1. - frtilgi)
    vern    = 0
    wal     = 1000. * rootdm * wci
    wapl    = 0
    waps    = 0
    was     = 0
    wetstor = 0

    # From parameters_plant.f90
    rehardredend    = 91.

    # From parameters_site.f90
    #    Simulation period and time step
    delt            = 1
    #    Atmospheric conditions
    co2a            = 350   
    #    Soil
    drate           = 50
    #    Soil - winter parameters
    lambdaice       = 1.9354e+005
    latentheat      = 335000.
    poolinfillimit  = 0.2
    rhowater        = 1000.
    #    Management: Irrigation Fraction
    irrigf          = 0
    #    Mathematical constants
    freq            = 2.*pi / 365.
    kmin            = 4.
    ampl            = 0.625
    bias            = kmin + ampl

    # Additional constants from translation
    # Needed to be defined once to allow variable to transit across loops
    freezepl = 0
    infil = 0
    packmelt = 0
    pooldrain = 0 
    poolinfil = 0 
    psnow = 0 
    refreeze = 0
    snowmelt = 0
    thawps = 0
    wremain = 0
    pwater = 0
    pinfil = 0
    runon = 0
    poolvolremain = 0
    eta = 0
    melt = 0
    wavail = 0
    tsurf = 0

    # Storage for PLOT output
    # E
    evp = []
    # T
    trn = []
    # Phenology
    pheno = 0
    ph = []
    # Tiller tilv
    tl_v = []
    # Non harvestable part of tillers (tilg1)
    tl_n = []
    # Harvestable part of tillers (tilg2)
    tl_h = []
    # Green Vegetation Biomass
    gv_b = []
    # Dead Vegetation Biomass
    dv_b = []
    # Harvested Biomass
    hv_Biomass = 0
    hv_b = []
    # LAI
    laim = []
    # NELLVG
    nell_vg = []

    # Compute the number of days of the run
    if( startyear == endyear):
        ndays = enddoy-startdoy
    else:
        ndays = (365 - startdoy) + enddoy + (endyear - startyear -1 ) * 365

    # Compute the starting date of the temporal index
    init_y_idx = [ i for i, y in enumerate(yeari) if y == startyear ]
    init_d_idx = [ i for i, d in enumerate(doyi) if d == startdoy ]
    # Intersect both indices list and get unique position
    init_idx = [ value for value in init_y_idx if value in init_d_idx ]
    if( startyear == endyear):
        stop_idx = ndays
    else:
        stop_idx = len(yeari) - ndays

    # Process the daily iteration from startdoy to enddoy across years
    for day in range(init_idx[0], stop_idx, 1):
        year    = yeari[day]    # day of the year (d)
        doy     = doyi[day]     # day of the year (d)
        rain    = raini[day]    # precipitation (mm d-1)
        gr      = gri[day]      # irradiation (mj m-2 d-1)
        tmmn    = tmmni[day]    # minimum (or average) temperature (degrees celsius)
        tmmx    = tmmxi[day]    # maximum (or average) temperature (degrees celsius)
        vp      = vpi[day]      # vapour pressure (kpa)
        wn      = wni[day]      # mean wind speed (m s-1)
        rsevap  = rsevapi[day]  # RS evaporation [optional, value is = 0.0 if not used]
        rstran  = rstrani[day]  # RS transpiration [optional, value is = 0.0 if not used]
        rslai   = rslaii[day]   # RS LAI [optional, value is = 0.0 if not used]
        rsgcut  = rsgcuti[day]  # RS grass cut event (0/1)

        # If RS data loaded, overwrite internal values with RS values
        if(rsevap != 0.0):
            evap = rsevap
        if(rstran != 0.0):
            tran = rstran
        if(rslai != 0.0):
            lai = rslai

        ############################################################################################
        # Compute equations
        ############################################################################################
        davtmp = (tmmn + tmmx) / 2.0
        dtr    = gr * exp( -ksnow * drystor)
        par    = 0.5 * 4.56 * dtr

        wcl = soilwatercontent(fdepth,rootd,wal)
        tsurf, fperm, frate = physics(davtmp, fdepth, rootd, sdepth, was, gamma, wcl, delt, wcfc, tsurf, lambdasoil, rhowater, latentheat)
	    # Microclimate
        pwater, psnow, snowmelt, wmaxstore, refreeze, staywet, wavail, wremain, wsupply, density, packmelt, rnintc, pinfil, runon, freezepl, thawps, pooldrain = rainsnowsurfacepool(davtmp, trainsnow, rain, bias, ampl, freq, doy, tmeltfreeze, melt, swret, swrf, drystor, snowmelt, wetstor, refreeze, pwater, wavail, sdepth, delt, rhopack, lai, fdepth, poolinfillimit, pinfil, wpoolmax, wapl, waps, runon, poolvolremain, poolinfil, frate, tsurf, lambdaice, rhowater, latentheat, eta, freezepl, infil, packmelt, pooldrain, thawps, wremain)
        if (waps == 0.): 
            permgas = 1.
        else: 
            permgas = 0.
        # End Microclimate
        # Day length
        dayl = ddayl(doy, lat)
        # Potential E and T
        pevap, ptran = penman(dtr, davtmp, wn, vp, lai, rnintc)
        # Light Interception
        parav, parint, dtrint = light(dayl, par, k, lai, dtr)
        # Actual E and T, evaporative fraction
        evap, tran, tranrf = evaptrtrf(fdepth, rootd, wal, wcad, pevap, wcfc, wcwp, ptran, tranco, wcst, wcwet, delt)
        # If we have RS data for both then replace them
        if(rsevap != 0 and rstran != 0):
            evap = rsevap               # (mm d-1)
            tran = rstran               # (mm d-1)
        # Root processes
        wcl, rrootd, explor = rootdg(fdepth, rootd, wal, rootdm, wcwp, rrdmax, delt, wcfc)
        # Soil Processes
        freezel, thaws, drain, runoff, irrig = frdrunir(wcfc, rootd, fdepth, wcst, infil, pooldrain, wal, delt, evap, tran, frate, was, drate, irrigf)
        # Root Oxygenation level
        fo2 = o2status(o2, rootd, fgas)
        # Plant processes
        #   Harvest if rsgcut != 0
        harv, noharv, fractv, harvla, harvlv, harvph, harvre, harvst, gstub, harvtilg2 = harvest(rsgcut, tilv, tilg2, claiv, lai, delt, clv, phen, hagere, cres, cst)
        #   Compute biomass
        cresmx, resnor = biomass(cocresmx, clv, cres, cst)
        #   Phenology
        gphen, dphen, phenrf, daylge = phenology(davtmp, daylp, dayl, daylb, phen, delt, phencr, dlmxge)
        print("Phenrf = %.2f, phen = %.2f, phencr = %.2f" % (phenrf, phen, phencr))
        #   Leaf growth
        efftmp, lerv, lerg, slamin, slanew = foliage1(tbase, davtmp, daylge, slamax, fslamin, resnor)
        #   LUE
        luemxq = lueco2tm(davtmp, rubisc, co2a, kluetilg, fractv, k, parav)
        #   Respiration Hardening Sink
        rehardredstart, resphardsi = hardeningsink(rehardredend, rehardredday, doy, tsurf, thardmx, lt50, lt50mn, hparam, clv, kresphard, resnor)
        #   Growth function
        nellvg, resmob, resphard, gres, glv, gst, grt, respgsh, respgrt = growth(parint, tranrf, luemxq, noharv, cres, tcres, davtmp, resphardsi, cresmx, delt, tilg2, cst, cstavm, simax1t, phenrf, nellvm, lerv, tilv, lfwidv, lerg, lfwidg, shape, slanew, yg, daylge)
        #   Plant Respiration
        rplantaer = plantrespiration(fo2, fo2mx, respgrt, respgsh, resphard)
        #   Senescence function
        dtanaer, hardrate, dehardrate, tv2til, dlai, dlv, dstub, dtilv, drt = senescence(permgas, tanaer, delt, ldt50a, ldt50b, lt50, krdranaer, krsr3h, tsurf, dparam, lt50mx, tsurfdiff, ratedmx, resphard, clv, kresphard, lai,laicr, rdrsco, rdrsmx, rdrtmin, rdrtem, noharv, cstub, rdrstub, tilv, crt, rdroot)
        #   Leaf Decline
        glai, gtilv, tilvg1, tilg1g2 = foliage2(slanew, glv, tsurf, tbase, phy, noharv, tranrf, daylge, fractv, phenrf, laitil, laieft, lai, fsmax, resnor, tilv, davtmp, toptge, tv2til, rgenmx, vern, dayl, daylg1g2, tilg1, rgrtg1g2)
        # Soil 2
        o2mx, o2out, o2in = o2fluxes(fo2mx, rootd, fgas, rplantaer, krtotaer, permgas, o2, delt)
        
        ############################################################################################
        # Update State variables
        ############################################################################################
        # Leaves
        clv     = clv + glv - dlv - harvlv
        # Dead Leaves
        clvd    = clvd + dlv
        # Harvested Biomass ((leaves + stems + stubble) / 0.45) + (reserves / 0.40)
        yielD   = (harvlv + harvst * hagere) / 0.45 + harvre / 0.40
        if (yielD > 0):
            yield_last = yielD
        # Respiration Balance (biomass?)
        cres    = cres + gres - resmob - harvre
        # Root Biomass
        crt     = crt + grt - drt
        # Stem Biomass
        cst     = cst + gst - harvst
        # Stubs Biomass
        cstub   = cstub + gstub - dstub
        # Transiting below zero water balance
        drystor = drystor + refreeze + psnow - snowmelt
        # Freezing Depth
        fdepth  = fdepth + frate
        # LAI generation (RS or computation)
        if(rslai != 0):
            lai = rslai
        else:
            lai = lai + glai - dlai - harvla
            #print("LAI = %.2f + %.2f - %.2f - %.2f" % (lai,glai,dlai,harvla))
        # Sensitivity to frost (Lethal Temperature for 50% of plants)
        lt50    = lt50 + dehardrate - hardrate
        # Oxygen balance
        o2      = o2 + o2in - o2out
        # Phenology
        phen    = min(1., phen + gphen - dphen - harvph)
        # Rooting depth
        rootd   = rootd + rrootd
        # Snow depth
        sdepth  = sdepth + (psnow / rhonewsnow) - packmelt
        # Anaerobic Transpiration
        tanaer  = tanaer + dtanaer
        # Tillers development
        #    Non-Harvestable part of tillers
        tilg1   = tilg1 + tilvg1 - tilg1g2
        #    Harvestable part of tillers
        tilg2   = tilg2 + tilg1g2 - harvtilg2
        #    Total Tillers
        tilv    = tilv + gtilv - tilvg1 - dtilv
        # Vernalisation
        if (davtmp < tvern):
            vern = 1
        # Water fluxes
        # Soil water amount: liquid (mm)
        wal     = wal + thaws - freezel + pooldrain + infil + explor + irrig - drain - runoff - evap - tran
        # Pool water amount: liquid (mm)
        wapl    = wapl + thawps - freezepl + poolinfil - pooldrain
        # Pool water amount: solid (=ice) (mm)
        waps    = waps - thawps + freezepl
        # Soil water amount: solid (=ice) (mm) 
        was     = was - thaws + freezel
        # Liquid water in snow (mm)
        wetstor = wetstor + wremain
        
        ############################################################################################
        # Append to PLOT lists
        ############################################################################################
        # NELLVG
        nell_vg.append(nellvg)
        #LAI
        laim.append(lai)
        # Phenology
        #pheno += phen
        ph.append(phen)
        # E
        evp.append(evap)
        # T
        trn.append(tran)
        # Tillers
        tl_v.append(tilv)
        #    Non-Harvestable part of tillers
        tl_n.append(tilg1)
        #    Harvestable part of tillers
        tl_h.append(tilg2)
        # Leaves + Stems (Green Vegetation Biomass)
        gv_b.append(clv+cst)
        # Dead [Leaves] (Dry Vegetation Biomass)
        dv_b.append(clvd)
        # Harvested cumulative Biomass
        hv_Biomass += yielD
        hv_b.append(hv_Biomass)

    return(nell_vg, laim, ph, evp, trn, tl_v, tl_n, tl_h, gv_b, dv_b, hv_b)
