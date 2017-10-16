import numpy as np

def rd_cdo(pth, mthd = 'bin_avg1day'):

    import csv
    from datetime import datetime
    # import pdb
    # import matplotlib.dates as mdates

    # Load data
    print 'Loading Juneau weather station data...'
    csvData = []
    csvHeader = []
    csvFileObj = open(pth)
    readerObj = csv.reader(csvFileObj)
    for row in readerObj:
        # pdb.set_trace()
        if row[0][:7] == 'STATION':
            csvHeader.append(row)    # read header
            continue
        csvData.append(row)
    csvFileObj.close()

    csvData = np.asarray(csvData)
    csvHeader = csvHeader[0]

    for i in range(np.size(csvHeader)):
        if csvHeader[i]=='STATION':
            stn = csvData[0, 0]
        elif csvHeader[i]=='STATION_NAME':
            stn_name = csvData[0, 1]
        elif csvHeader[i]=='ELEVATION':
            elv = np.mean(csvData[:, i].astype(np.float))
        elif csvHeader[i]=='LATITUDE':
            lat = np.mean(csvData[:, i].astype(np.float))
        elif csvHeader[i]=='LONGITUDE':
            lon = np.mean(csvData[:, i].astype(np.float))
        elif csvHeader[i]=='DATE':
            pyt = np.array([datetime.strptime(csvData[j, i], '%Y%m%d') for j in range(csvData[:, i].size)])
        elif csvHeader[i]=='PRCP':
            prcp = csvData[:, i].astype(np.float)  # [mm]
        elif csvHeader[i]=='SNWD':
            snwd = csvData[:, i].astype(np.float)  # [mm]
        elif csvHeader[i]=='SNOW':
            snow = csvData[:, i].astype(np.float)  # [mm]
        elif csvHeader[i]=='ACMH':
            acmh = csvData[:, i].astype(np.float)  # [percent]
        elif csvHeader[i]=='ACSH':
            acsh = csvData[:, i].astype(np.float)  # [percent]
        elif csvHeader[i]=='TSUN':
            tsun = csvData[:, i].astype(np.float)  # [minutes]
        elif csvHeader[i]=='TMAX':
            tmax = csvData[:, i].astype(np.float)  # [degC]
        elif csvHeader[i]=='TMIN':
            tmin = csvData[:, i].astype(np.float)  # [degC]
        elif csvHeader[i]=='AWND':
            awnd = csvData[:, i].astype(np.float)  # [ms-1]
        elif csvHeader[i]=='WDF1':
            wdf1 = csvData[:, i].astype(np.float)  # [deg]
        elif csvHeader[i]=='WDF2':
            wdf2 = csvData[:, i].astype(np.float)  # [deg]
        elif csvHeader[i]=='WDF5':
            wdf5 = csvData[:, i].astype(np.float)  # [deg]
        elif csvHeader[i]=='WSF1':
            wsf1 = csvData[:, i].astype(np.float)  # [ms-1]
        elif csvHeader[i]=='WSF2':
            wsf2 = csvData[:, i].astype(np.float)  # [ms-1]
        elif csvHeader[i]=='WSF5':
            wsf5 = csvData[:, i].astype(np.float)  # [ms-1]

    # Mask out invalid data
    print 'Setting invalid data to NaN...'
    prcp[prcp==-9999]  = np.NaN
    snwd[snwd==-9999]  = np.NaN
    snow[snow==-9999]  = np.NaN
    acmh[acmh==-9999]  = np.NaN
    acsh[acsh==-9999]  = np.NaN
    tsun[tsun==-9999]  = np.NaN
    tmax[tmax==-9999]  = np.NaN
    tmin[tmin==-9999]  = np.NaN
    awnd[awnd==-9999]  = np.NaN
    wdf1[wdf1==-9999]  = np.NaN
    wdf2[wdf2==-9999]  = np.NaN
    wdf5[wdf5==-9999]  = np.NaN
    wsf1[wsf1==-9999]  = np.NaN
    wsf2[wsf2==-9999]  = np.NaN
    wsf5[wsf5==-9999]  = np.NaN
       
    # pdb.set_trace()
    # Mask out invalid data
    print 'Masking invalid data...'
    prcp = np.ma.masked_invalid(prcp)
    snwd = np.ma.masked_invalid(snwd)
    snow = np.ma.masked_invalid(snow)
    acmh = np.ma.masked_invalid(acmh)
    acsh = np.ma.masked_invalid(acsh)
    tsun = np.ma.masked_invalid(tsun)
    tmax = np.ma.masked_invalid(tmax)
    tmin = np.ma.masked_invalid(tmin)
    awnd = np.ma.masked_invalid(awnd)
    wdf1 = np.ma.masked_invalid(wdf1)
    wdf2 = np.ma.masked_invalid(wdf2)
    wdf5 = np.ma.masked_invalid(wdf5)
    wsf1 = np.ma.masked_invalid(wsf1)
    wsf2 = np.ma.masked_invalid(wsf2)
    wsf5 = np.ma.masked_invalid(wsf5)
    
    data = {'stn': stn,
            'stn_name': stn_name,
            'pyt': pyt,
            'prcp': prcp,
            'snwd': snwd,
            'snow': snow,
            'acmh': acmh,
            'acsh': acsh,
            'tsun': tsun,
            'tmax': tmax,
            'tmin': tmin,
            'awnd': awnd,
            'wdf1': wdf1,
            'wdf2': wdf2,
            'wdf5': wdf5,
            'wsf1': wsf1,
            'wsf2': wsf2,
            'wsf5': wsf5
           }

    return data


# The five core values are:
# PRCP = Precipitation (mm or inches as per user preference, inches to hundredths on Daily Form pdf file) SNOW = Snowfall (mm or inches as per user preference, inches to tenths on Daily Form pdf file)
# SNWD = Snow depth (mm or inches as per user preference, inches on Daily Form pdf file)
# TMAX = Maximum temperature (Fahrenheit or Celsius as per user preference, Fahrenheit to tenths on Daily Form pdf file
# TMIN = Minimum temperature (Fahrenheit or Celsius as per user preference, Fahrenheit to tenths on Daily Form pdf file
# The other values are:
# ACMC = Average cloudiness midnight to midnight from 30-second ceilometer data (percent) ACMH = Average cloudiness midnight to midnight from manual observations (percent) ACSC = Average cloudiness sunrise to sunset from 30-second ceilometer data (percent) ACSH = Average cloudiness sunrise to sunset from manual observations (percent)
# AWND = Average daily wind speed (meters per second or miles per hour as per user preference) DAEV = Number of days included in the multiday evaporation total (MDEV)
# DAPR = Number of days included in the multiday precipitation total (MDPR)
# DASF = Number of days included in the multiday snowfall total (MDSF)
# DATN = Number of days included in the multiday minimum temperature (MDTN) DATX = Number of days included in the multiday maximum temperature (MDTX) DAWM = Number of days included in the multiday wind movement (MDWM)
# DWPR = Number of days with non-zero precipitation included in multiday precipitation total (MDPR) EVAP = Evaporation of water from evaporation pan (mm or inches as per user preference, or hundredths of inches on Daily Form pdf file)
# FMTM = Time of fastest mile or fastest 1-minute wind (hours and minutes, i.e., HHMM)
# FRGB = Base of frozen ground layer (cm or inches as per user preference)
# FRGT = Top of frozen ground layer (cm or inches as per user preference)
# FRTH = Thickness of frozen ground layer (cm or inches as per user preference)
# GAHT = Difference between river and gauge height (cm or inches as per user preference)
# MDEV = Multiday evaporation total (mm or inches as per user preference; use with DAEV)
# MDPR = Multiday precipitation total (mm or inches as per user preference; use with DAPR and DWPR, if available)
# MDSF = Multiday snowfall total (mm or inches as per user preference)
# MDTN = Multiday minimum temperature (Fahrenheit or Celsius as per user preference ; use with DATN) MDTX = Multiday maximum temperature (Fahrenheit or Celsius as per user preference ; use with DATX) MDWM = Multiday wind movement (miles or km as per user preference)
# MNPN = Daily minimum temperature of water in an evaporation pan (Fahrenheit or Celsius as per user preference)
# MXPN = Daily maximum temperature of water in an evaporation pan (Fahrenheit or Celsius as per user preference)
# 
# PGTM = Peak gust time (hours and minutes, i.e., HHMM)
# PSUN = Daily percent of possible sunshine (percent)
# SN*# = Minimum soil temperature where * corresponds to a code
# for ground cover and # corresponds to a code for soil depth (Fahrenheit or Celsius as per user preference)
# Ground cover codes include the following: 0 = unknown
# 1 = grass
# 2 = fallow
# 3 = bare ground 4 = brome grass 5 = sod
# 6 = straw mulch 7 = grass muck 8 = bare muck
# Depth codes include the following: 1 = 5 cm
# 2 = 10 cm 3 = 20 cm 4 = 50 cm 5 = 100 cm 6 = 150 cm 7 = 180 cm
# SX*# = Maximum soil temperature where * corresponds to a code for ground
# cover and # corresponds to a code for soil depth. See SN*# for depth codes. (Fahrenheit or
# Celsius as per user preference)
# THIC = Thickness of ice on water (inches or mm as per user preference)
# TOBS = Temperature at the time of observation (Fahrenheit or Celsius as per user preference) TSUN = Daily total sunshine (minutes)
# WDF1 = Direction of fastest 1-minute wind (degrees)
# WDF2 = Direction of fastest 2-minute wind (degrees)
# WDF5 = Direction of fastest 5-second wind (degrees) WDFG = Direction of peak wind gust (degrees)
# WDFI = Direction of highest instantaneous wind (degrees) WDFM = Fastest mile wind direction (degrees)
# WDMV = 24-hour wind movement (km or miles as per user preference, miles on Daily Form pdf file) WESD = Water equivalent of snow on the ground (inches or mm as per user preference)
# WESF = Water equivalent of snowfall (inches or mm as per user preference)
# WSF1 = Fastest 1-minute wind speed (miles per hour or meters per second as per user preference) WSF2 = Fastest 2-minute wind speed (miles per hour or meters per second as per user preference) WSF5 = Fastest 5-second wind speed (miles per hour or meters per second as per user preference) WSFG = Peak guest wind speed (miles per hour or meters per second as per user preference)
# WSFI = Highest instantaneous wind speed (miles per hour or meters per second as per user preference) WSFM = Fastest mile wind speed (miles per hour or meters per second as per user preference)
# WT** = Weather Type where ** has one of the following values:
# 01 = Fog, ice fog, or freezing fog (may include heavy fog)
# 02 = Heavy fog or heaving freezing fog (not always distinguished from fog)
# 03 = Thunder
# 04 = Ice pellets, sleet, snow pellets, or small hail
# 05 = Hail (may include small hail)
# 06 = Glaze or rime
# 07 = Dust, volcanic ash, blowing dust, blowing sand, or blowing obstruction 08 = Smoke or haze
# 09 = Blowing or drifting snow
# 10 = Tornado, waterspout, or funnel cloud
# 11 = High or damaging winds
# 12 = Blowing spray
# 13 = Mist
# 14 = Drizzle
# 15 = Freezing drizzle
# 16 = Rain (may include freezing rain, drizzle, and freezing drizzle) 17 = Freezing rain
# 18 = Snow, snow pellets, snow grains, or ice crystals
# 19 = Unknown source of precipitation
# 21 = Ground fog
# 22 = Ice fog or freezing fog


