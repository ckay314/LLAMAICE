import numpy as np
import datetime
import matplotlib.pyplot as plt
import sys
import matplotlib.dates as mdates


obsDataPath = './ISdata/'

shadeICE = True
plotCats = True
plotHSS = True
addLabels = True

# ICE region colors
cols = ['r', 'orange', 'b', 'b', 'orange']     

# boundary catalog colors
catCols = {'ATSB':'Cyan', 'CDAW':'#FCB105', 'DONKI':'Maroon', 'DREAMS':'Magenta', 'GMU1':'#F0E68C', 'GMU2':'Gold', 'H4C':'#6915CE', 'Lepping':'#F39C12', 'RC':'Blue', 'Wind':'Red'}

# HSS catalog styles [plot offset, color]
catStyles = {'G':[-0.1, 'r'] , 'V':[0.,'b'], 'X':[0.,'orange'], 'D':[.1,'maroon']}

# Observation aesthetics
obsLW = 1.5 # linewidth 
ACEcol = 'k'
Windcol = 'gray'



# The other CME catalogs
otherData = np.genfromtxt('revICEless.dat', dtype=str)
# LLAMAICE results which are separate bc more bounds
ICEdata   = np.genfromtxt('LLAMAICE0.2.csv', delimiter=',', dtype=str)


CMEchoice = sys.argv[1]
ICEcase = True
failCase = False

totalCMEs = otherData[-1,0]

# ============================== ID whether ICE case or generic ==============================
# =====  check if date or CME number =====
# assume given a date if has - in it
if '-' in CMEchoice: 
    if CMEchoice in ICEdata[:,2]:
        myID = np.where(ICEdata[:,2] == CMEchoice)[0][0]
        myoIDs = np.where(otherData[:,1] == CMEchoice)[0]
        print (CMEchoice + ' corresponds to LLAMAICE case '+ICEdata[myID,0])
        # Get the appropriate data files
        myYr = ICEdata[myID,2][0:4]
        myDate = ICEdata[myID,2]
        myDT   = datetime.datetime.strptime(myDate, "%Y-%m-%dT%H:%M" )
    else:
        print ('Not an LLAMAICE CME. Plotting anyway...')
        try:
            # Get the appropriate data files
            myYr = CMEchoice[0:4]
            myDate = CMEchoice
            myDT   = datetime.datetime.strptime(CMEchoice, "%Y-%m-%dT%H:%M" )
            ICEcase = False
        except:
            print ('Date format error. Make sure in format YYYY-MM-DDTHH:MM ')
            sys.exit()
            
# assume given ICE number otherwise 
else:
    if int(CMEchoice) <= int(totalCMEs):
        myID = np.where(ICEdata[:,0].astype(int) == int(CMEchoice))[0][0]
        myoIDs = np.where(otherData[:,0] == CMEchoice)[0]
        # Get the appropriate data files
        dateStr = ICEdata[myID,2] # sheath time
        if dateStr == '-':
            dateStr = ICEdata[myID,4] # FR time
        if dateStr == '-': # fail/supsect case
            failCase = True
            print('No LLAMAICE bounds only other catalog(s).')
            dateStr = otherData[myoIDs[0],1]
        myYr = dateStr[0:4]
        print ('Results for CME '+dateStr)
        myDT   = datetime.datetime.strptime(dateStr, "%Y-%m-%dT%H:%M" )
    else:
        print ('Incorrect CME id. LLAMAICE currently contains '+totalCMEs+' CMEs.')
        sys.exit()



# ============================== Observations Processing ==============================
# ===== pull in obs data and format ===== 
hasACE = True
if int(myYr) < 1997:
    hasACE = False
if hasACE:
    aceDataFull = np.genfromtxt(obsDataPath+'ace5/ace_5min_'+myYr+'.dat', dtype=str)
windDataFull = np.genfromtxt(obsDataPath+'wind5/wind_MFISWE_5min_'+myYr+'.dat', dtype=str)

# Downselect to # days before and # after
startPlotDT = myDT - datetime.timedelta(days=2)
startStr    = startPlotDT.strftime("%Y-%m-%dT00:00")  
endPlotDT   = myDT + datetime.timedelta(days=4)
endStr    = endPlotDT.strftime("%Y-%m-%dT23:55")  

# Check if we are going over a year boundary
if startPlotDT.year != myDT.year:
    if int(myYr) >= 1998:
        bonusACE = np.genfromtxt(obsDataPath+'ace5/ace_5min_'+str(int(myYr)-1)+'.dat', dtype=str)
        aceDataFull = np.concatenate((bonusACE, aceDataFull), axis=0)
    if int(myYr) >= 1996:
        bonusWind = np.genfromtxt(obsDataPath+'wind5/wind_MFISWE_5min_'+str(int(myYr)-1)+'.dat', dtype=str)
        windDataFull = np.concatenate((windDataFull, bonusWind), axis=0)
elif endPlotDT.year != myDT.year:
    if int(myYr) <= 2022:
        bonusACE = np.genfromtxt(obsDataPath+'ace5/ace_5min_'+str(int(myYr)+1)+'.dat', dtype=str)
        aceDataFull = np.concatenate((bonusACE, aceDataFull), axis=0)
        bonusWind = np.genfromtxt(obsDataPath+'wind5/wind_MFISWE_5min_'+str(int(myYr)+1)+'.dat', dtype=str)
        windDataFull = np.concatenate(( windDataFull, bonusWind), axis=0)

# Down select to nice region
if hasACE:
    try:
        startidx = np.where(aceDataFull[:,0] == startStr)[0][0]
    except:
        startidx = 0
    try:
        endidx  = np.where(aceDataFull[:,0] == endStr)[0][0]
    except:
        endidx = len(aceDataFull[:,0])
    aceData = aceDataFull[startidx:endidx+1, :]

try:    
    startidx = np.where(windDataFull[:,0] == startStr)[0][0]
except:
    startidx = 0 
try:
    endidx  = np.where(windDataFull[:,0] == endStr)[0][0]
except:
    endidx = len(windDataFull[:,0])
windData = windDataFull[startidx:endidx+1, :]

# make a date time array
if hasACE:
    aceTimes = []
    for i in range(len(aceData)):
        newDT = datetime.datetime.strptime(aceData[i,0], "%Y-%m-%dT%H:%M" )
        aceTimes.append(newDT)
    aceTimes = np.array(aceTimes)   

windTimes = []
for i in range(len(windData)):
    newDT = datetime.datetime.strptime(windData[i,0], "%Y-%m-%dT%H:%M" )
    windTimes.append(newDT)
windTimes = np.array(windTimes)   

 
# ===== Caclulate things and clean up data =====
if hasACE:
    BidxA = np.where((aceData[:,1] != '-9999') & (aceData[:,1].astype(float) > -999.))[0]
    BxACE, ByACE, BzACE =  aceData[:,1].astype(float),  aceData[:,2].astype(float),  aceData[:,3].astype(float)
    BtotACE = np.sqrt(BxACE**2 + ByACE**2 + BzACE**2)
    BthetaA = np.arctan2(BzACE, np.sqrt(BxACE**2+ByACE**2))*180/np.pi
    BphiA = np.arctan2(-ByACE, -BxACE)*180/np.pi   
    BphiA[BphiA<0] =  BphiA[BphiA<0] + 360.
    vsACE = aceData[:,4].astype(float)
    TsACE = aceData[:,5].astype(float)
    nsACE = aceData[:,6].astype(float)
    vidxA = np.where((vsACE < 10000) & (vsACE > 0))[0]
    TidxA = np.where((TsACE < 5e6) & (TsACE > 1000))[0]
    nidxA = np.where((nsACE < 200) & (nsACE > 0))[0]
    betaIdxA = []
    for idx in BidxA:
        if idx in TidxA:
            betaIdxA.append(idx)    
    # calc plasma beta
    betaACE = (4.16e-5 * TsACE[betaIdxA] + 5.34) * nsACE[betaIdxA]/BtotACE[betaIdxA]**2 # omni version of calc w/Tp, np, dunno where 5.34 comes from
    
BidxW = np.where((windData[:,1] != '-9999'))[0]
BxWind, ByWind, BzWind =  windData[:,1].astype(float),  windData[:,2].astype(float),  windData[:,3].astype(float)
BtotWind = np.sqrt(BxWind**2 + ByWind**2 + BzWind**2)
BidxW = np.where((windData[:,1] != '-9999') & (BtotWind < 100))[0]

BthetaW = np.arctan2(BzWind, np.sqrt(BxWind**2+ByWind**2))*180/np.pi
BphiW = np.arctan2(-ByWind, -BxWind)*180/np.pi        
BphiW[BphiW<0] =  BphiW[BphiW<0] + 360.
vsWind = windData[:,4].astype(float)
TsWind = windData[:,5].astype(float)
nsWind = windData[:,6].astype(float)
vidxW = np.where((vsWind < 10000) & (vsWind > 0))[0]
TidxW = np.where(TsWind < 5e6)[0]
nidxW = np.where((nsWind < 200) & (nsWind  > 0))[0]
betaIdxW = []
for idx in BidxW:
    if idx in TidxW:
        betaIdxW.append(idx)
# calc plasma beta
betaWind = (4.16e-5 * TsWind[betaIdxW] + 5.34) * nsWind[betaIdxW]/BtotWind[betaIdxW]**2 # omni version of calc w/Tp, np, dunno where 5.34 comes from



# ============================== Read in ICE boundaries ==============================
toPlot = []
ICEbounds = None
bonusLines = []
if ICEcase:
    # Check all lines with this ID number
    for idx in myoIDs:
        myVals = otherData[idx,3], otherData[idx,4], otherData[idx,5], otherData[idx,6] 
        # Get the other catalog boundaries
        if (myVals not in toPlot) & (otherData[idx,3] != 'LLAMAICE'):
            toPlot.append(myVals)
        # Ignore duplicates (though all rm'ed now?)
        else:
           print ('Duplicate for '+myVals[0])
    # Save ICE separate
    ICEbounds = 'LLAMAICE', ICEdata[myID, 2], ICEdata[myID, 3], ICEdata[myID, 4], ICEdata[myID, 5], ICEdata[myID, 6],  ICEdata[myID, 7]

# ============================== Plotting action ==============================
def mouse_event(event):
    clickstr = mdates.num2date(event.xdata).strftime("%Y-%m-%dT%H:%M")
    clickDT = datetime.datetime.strptime(clickstr, "%Y-%m-%dT%H:%M" )
    doy = (clickDT - datetime.datetime(clickDT.year, 1, 1)).total_seconds()/24/3600.
    print('Date: {} and y-val: {}, doy: {}'.format(clickstr, event.ydata, '{:6.3f}'.format(doy+1.) ))

    
# ============================== Set up and make figure ============================== 
fig, ax0 = plt.subplots(11, 1, sharex=True, figsize=(10,7.5), gridspec_kw = {'height_ratios':[0.3,1,1,1,1,1,1,1,1,1,1]}) 
axes = [ax0[1], ax0[2], ax0[3],ax0[4],ax0[5],ax0[6],ax0[7],ax0[8],ax0[9],ax0[10],ax0[0]]
cid = fig.canvas.mpl_connect('button_press_event', mouse_event)


# Plot ACE if exists
if hasACE:
    axes[0].plot(aceTimes[BidxA], BtotACE[BidxA], ACEcol, lw=obsLW)
    axes[1].plot(aceTimes[BidxA], BxACE[BidxA], ACEcol, lw=obsLW)
    axes[2].plot(aceTimes[BidxA], ByACE[BidxA], ACEcol, lw=obsLW)
    axes[3].plot(aceTimes[BidxA], BzACE[BidxA], ACEcol, lw=obsLW)
    axes[4].plot(aceTimes[BidxA], BthetaA[BidxA], ACEcol, lw=obsLW)
    axes[5].plot(aceTimes[BidxA], BphiA[BidxA], ACEcol, lw=obsLW)
    axes[6].plot(aceTimes[vidxA], vsACE[vidxA], ACEcol, lw=obsLW)
    axes[7].plot(aceTimes[TidxA], TsACE[TidxA], ACEcol, lw=obsLW)
    axes[8].plot(aceTimes[nidxA], nsACE[nidxA], ACEcol, lw=obsLW)
    axes[9].plot(aceTimes[betaIdxA], betaACE, ACEcol, lw=obsLW)

# Assume has Wind, plot it
axes[0].plot(windTimes[BidxW], BtotWind[BidxW], Windcol, lw=obsLW)
axes[1].plot(windTimes[BidxW], BxWind[BidxW], Windcol, lw=obsLW)
axes[2].plot(windTimes[BidxW], ByWind[BidxW], Windcol, lw=obsLW)
axes[3].plot(windTimes[BidxW], BzWind[BidxW], Windcol, lw=obsLW)
axes[4].plot(windTimes[BidxW], BthetaW[BidxW], Windcol, lw=obsLW)
axes[5].plot(windTimes[BidxW], BphiW[BidxW], Windcol, lw=obsLW)
axes[6].plot(windTimes[vidxW], vsWind[vidxW], Windcol, lw=obsLW)
axes[7].plot(windTimes[TidxW], TsWind[TidxW], Windcol, lw=obsLW)
axes[7].set_yscale('log')
axes[8].plot(windTimes[nidxW], nsWind[nidxW], Windcol, lw=obsLW)
axes[9].plot(windTimes[betaIdxW], betaWind, Windcol, lw=obsLW)
axes[9].set_yscale('log')
# add line at beta = 1
xl = axes[9].get_xlim()
axes[9].plot(xl, [1,1], 'k--', zorder=0 )
axes[9].set_xlim(xl)    
                   
# get boundaries at end of obs plotting
bounds = []
for i in range(10):
    bounds.append(axes[i].get_ylim())
 
 
 
# ============================== Polishing plot limits ============================== 
defMins = [0, -15,  -15, -15, -90, 0, 300, 1e4, 0, 1e-2]
defMaxs = [25, 15, 15,  15, 90, 360, 800, 1e6, 25, 3]
minFlags = [0, -300, -300, -300, -90, 0, 0, 100, 0, 1e-4]
maxFlags = [1000, 300, 300, 300, 90, 360, 4000, 5e6,150, 5e6]
# set time bound for ICE cases and reset y ranges based on subset of obs data
if ICEcase:
    plotbounds = np.genfromtxt('plotBounds.dat', dtype=str)
    bidx = np.where(plotbounds[:,0] == CMEchoice)[0]
    pstart = datetime.datetime.strptime(plotbounds[bidx[0],2], "%Y-%m-%dT%H:%M" )
    pend   = datetime.datetime.strptime(plotbounds[bidx[0],3], "%Y-%m-%dT%H:%M" )
    axes[-1].set_xlim([pstart, pend])
    
    # reset y lims based on current range ------------------------------------------------------------
    windIdx = np.where((windTimes >= pstart) & (windTimes <= pend))
    windMins = np.zeros(10)
    windMaxs = np.zeros(10)
    allWind = [BtotWind, BxWind, ByWind, BzWind, BthetaW, BphiW, vsWind, TsWind, nsWind, betaWind]
    for i in range(10):
        # try expect makes happy when one param has data removed
        flagIt = False
        try:
            windMins[i] = np.min(allWind[i][windIdx])
        except:
            try:
                windMins[i] = np.min(allWind[i])
            except:
                windMins[i] = defMins[i]
            flagIt = True
        try:
            windMaxs[i] = np.max(allWind[i][windIdx])
        except:
            try:
                windMaxs[i] = np.max(allWind[i])
            except:
                windMaxs[i] = defMaxs[i]
            flagIt = True
            
        if windMins[i] == -9999: 
            subvals = allWind[i][windIdx]
            try:
                windMins[i] = np.min(subvals[np.where(subvals != -9999)])
            except:
                windMins[i] = defMins[i]
        if windMaxs[i] == -9999:
            subvals = allWind[i][windIdx]
            try:
                windMaxs[i] = np.max(subvals[np.where(subvals != -9999)])
            except:
                windMaxs[i] = defMaxs[i]

    if hasACE:
        allACE = [BtotACE, BxACE, ByACE, BzACE, BthetaA, BphiA, vsACE, TsACE, nsACE, betaACE]
        aceIdx = np.where((aceTimes >= pstart) & (aceTimes <= pend))
        aceMins = np.zeros(10)
        aceMaxs = np.zeros(10)
        for i in range(10):
            try:
                aceMins[i] = np.min(allACE[i][aceIdx])
            except:
                aceMins[i] = 999999
            try:
                aceMaxs[i] = np.max(allACE[i][aceIdx])
            except:
                aceMaxs[i] = -999999
            if aceMins[i] == -999.9:
                goodIdx = np.where(allACE[i][aceIdx] != -999.9)
                aceMins[i] = np.min(allACE[i][aceIdx][goodIdx])
                if aceMins[i] == -999.9:
                    aceMins[i] = windMins[i]
            if aceMaxs[i] == -999.9:
                goodIdx = np.where(allACE[i][aceIdx] != -999.9)
                aceMaxs[i] = np.max(allACE[i][aceIdx][goodIdx])
                if aceMaxs[i] == -999.9:
                    aceMaxs[i] = windMaxs[i]
    else:
        aceMins = windMins
        aceMaxs = windMaxs
    
    allMins = np.zeros(10)
    allMaxs = np.zeros(10)
    for i in range(10):
        allMins[i] = np.min([aceMins[i], windMins[i]])
        allMaxs[i] = np.max([aceMaxs[i], windMaxs[i]])
        if (allMaxs[i] > 10 * np.min([aceMaxs[i], windMaxs[i]])) & (np.min([aceMaxs[i], windMaxs[i]]) >0) & hasACE:
            # try checking the 99th percentile of the other sat to see if rm bad data
            try:
                nearMax = [np.percentile(allWind[i][windIdx],99),np.percentile(allACE[i][aceIdx],99)]
                if np.max(nearMax) < 2*np.min(nearMax):
                    allMaxs[i] = np.max(nearMax) 
                else:
                    allMaxs[i] = np.min([aceMaxs[i], windMaxs[i]])
            except:
                allMaxs[i] = np.min([aceMaxs[i], windMaxs[i]])
        elif (aceMaxs[i] == windMaxs[i]) & (allMaxs[i] > defMaxs[i]):
            # assume we only have Wind and it is a little wonky
            my99 = np.percentile(allWind[i][windIdx], 99)
            if my99 < defMaxs[i]:
                allMaxs[i] = my99
        if allMaxs[i] > bounds[i][1]:
            allMaxs[i] = bounds[i][1]
        if (allMins[i] < 0) & (allMins[i] < 5 * np.max([aceMins[i], windMins[i]])):
            if np.max([aceMins[i], windMins[i]]) < -1:
                allMins[i] = np.max([aceMins[i], windMins[i]])
        if (allMins[i] > 0) & (allMins[i] < 5) & (i != 9):
            allMins[i] = 0
        
        if allMins[i] < minFlags[i]:
            otherMin = np.max([aceMins[i], windMins[i]])
            if (otherMin > minFlags[i]) & (otherMin < maxFlags[i]):
                allMins[i] = otherMin
            else:
                allMins[i] = minFlags[i]
                
        if allMaxs[i] > maxFlags[i]:
            otherMax = np.min([aceMaxs[i], windMaxs[i]])
            if (otherMax < maxFlags[i]) & (otherMax > minFlags[i]):
                allMaxs[i] = otherMax
            else:
                allMaxs[i] = maxFlags[i]
                        
        #print (i, allMins[i], aceMins[i], windMins[i], allMaxs[i], aceMaxs[i], windMaxs[i])
        if i in [0,5,7,8,9]: # min near zero
            axes[i].set_ylim(0.8*allMins[i], 1.1*allMaxs[i])
        elif i == 6:
            axes[i].set_ylim(allMins[i]-50, allMaxs[i]+50)
        else: # min is negative
            axes[i].set_ylim(1.1*allMins[i], 1.1*allMaxs[i])
        
     
    # ------- Collect in situ boundaries from nearby cases --------
    nBonus = 20
    moreBounds = []
    try:
        for i in range(nBonus):
            # Front side
            t1 = otherData[myoIDs[0]-i-1, 4]
            t2 = otherData[myoIDs[0]-i-1, 5]
            t3 = otherData[myoIDs[0]-i-1, 6]
            if t1 != 'None':
                t1 = datetime.datetime.strptime(t1, "%Y-%m-%dT%H:%M:%S" )
                if t1 > pstart:
                    moreBounds.append(t1)
            if t2 != 'None':
                t2 = datetime.datetime.strptime(t2, "%Y-%m-%dT%H:%M:%S" )
                if t2 > pstart:
                    moreBounds.append(t2)
            if t3 != 'None':
                t3 = datetime.datetime.strptime(t3, "%Y-%m-%dT%H:%M:%S" )
                if t3 > pstart:
                    moreBounds.append(t3)
    except:
        pass # assume exit bc hitting idx 0
    try:
        for i in range(nBonus):
            # Back side
            t4 = otherData[myoIDs[-1]+i+1, 4]
            t5 = otherData[myoIDs[-1]+i+1, 5]
            t6 = otherData[myoIDs[-1]+i+1, 6]
            if t4 != 'None':
                t4 = datetime.datetime.strptime(t4, "%Y-%m-%dT%H:%M:%S" )
                if t4 < pend:
                    moreBounds.append(t4)
            if t5 != 'None':
                t5 = datetime.datetime.strptime(t5, "%Y-%m-%dT%H:%M:%S" )
                if t5 < pend:
                    moreBounds.append(t5)
            if t6 != 'None':
                t6 = datetime.datetime.strptime(t6, "%Y-%m-%dT%H:%M:%S" )
                if t6 < pend:
                    moreBounds.append(t6)
    except:
        pass # assume exit bc hit end
        
    # add all the nearby in situ boundaries
    bounds = []
    for i in range(10):
        bounds.append(axes[i].get_ylim())
    for item in moreBounds:
        for i in range(10):
            axes[i].plot([item, item], bounds[i], 'gray', linestyle=':')
    for i in range(10):
        axes[i].set_ylim(bounds[i])

# ============================== Plot the boundaries for this event ============================== 
if plotCats:
    # get boundaries
    bounds = []
    for i in range(10):
        bounds.append(axes[i].get_ylim())
    printNames = []
    print('')
    print ('Catalogs including this event:')
    for item in toPlot:
        if item[0] not in printNames:
            printNames.append(item[0])
        if item[1] != 'None':
            for i in range(10):
                thisDT = datetime.datetime.strptime(item[1], "%Y-%m-%dT%H:%M:%S" )
                axes[i].plot([thisDT, thisDT], bounds[i], catCols[item[0]], linestyle='--')
        if item[2] != 'None':
            for i in range(10):
                thisDT = datetime.datetime.strptime(item[2], "%Y-%m-%dT%H:%M:%S" )
                axes[i].plot([thisDT, thisDT], bounds[i], catCols[item[0]], linestyle='--')
        if item[3] != 'None':
            for i in range(10):
                thisDT = datetime.datetime.strptime(item[3], "%Y-%m-%dT%H:%M:%S" )
                axes[i].plot([thisDT, thisDT], bounds[i], catCols[item[0]], linestyle='--')
        print (item[0].rjust(11)  + item[1].rjust(23) + item[2].rjust(23) + item[3].rjust(23))


# ==============================  Shade in LLAMAICE regions ==============================  
if ICEbounds:
    alp = 0.2
    if shadeICE:
        for i in range(10):
            if (ICEbounds[1] != '-') & (ICEbounds[3] != '-'):
                if ICEbounds[2] in ['-', '']:
                    ib1 = datetime.datetime.strptime(ICEbounds[1], "%Y-%m-%dT%H:%M" )
                    ib2 = datetime.datetime.strptime(ICEbounds[3], "%Y-%m-%dT%H:%M" )
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[0], zorder=0, alpha=alp )
                else:
                    ib1 = datetime.datetime.strptime(ICEbounds[1], "%Y-%m-%dT%H:%M" )
                    ib2 = datetime.datetime.strptime(ICEbounds[2], "%Y-%m-%dT%H:%M" )
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[0], zorder=0, alpha=alp )
                    ib1 = datetime.datetime.strptime(ICEbounds[2], "%Y-%m-%dT%H:%M" )
                    ib2 = datetime.datetime.strptime(ICEbounds[3], "%Y-%m-%dT%H:%M" )
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[1], zorder=0, alpha=alp )
            elif (ICEbounds[1] == '-') & (ICEbounds[2] not in ['-','']):   
                ib1 = datetime.datetime.strptime(ICEbounds[2], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(ICEbounds[3], "%Y-%m-%dT%H:%M" )
                axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[1], zorder=0, alpha=alp )         
            if (ICEbounds[3] != '-') & (ICEbounds[4] != '-'):
                ib1 = datetime.datetime.strptime(ICEbounds[3], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(ICEbounds[4], "%Y-%m-%dT%H:%M" )
                axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[2], zorder=0, alpha=alp )
            if (ICEbounds[4] != '-') & (ICEbounds[5] != '-'):
                ib1 = datetime.datetime.strptime(ICEbounds[4], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(ICEbounds[5], "%Y-%m-%dT%H:%M" )
                axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[2], zorder=0, alpha=alp/2 )
            if (ICEbounds[5] != '-') & (ICEbounds[6] != '-'):
                ib1 = datetime.datetime.strptime(ICEbounds[5], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(ICEbounds[6], "%Y-%m-%dT%H:%M" )
                axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[4], zorder=0, alpha=alp )
            if (ICEbounds[5] == '-') & (ICEbounds[4] != '-') & (ICEbounds[6] != '-'):
                ib1 = datetime.datetime.strptime(ICEbounds[4], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(ICEbounds[6], "%Y-%m-%dT%H:%M" )
                axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[4], zorder=0, alpha=alp)
    print (ICEbounds[0].rjust(11)  + ICEbounds[1].rjust(23) + ICEbounds[2].rjust(23) + ICEbounds[3].rjust(23)+ ICEbounds[4].rjust(23)+ ICEbounds[5].rjust(23)+ ICEbounds[6].rjust(23))
print('')            


# ============================== Shade in  adjacent event ============================== 
if ICEcase:
    alp =  0.07
    # shade in prev case
    if (myID !=0) & shadeICE:
        prevICE = 'LLAMAICE', ICEdata[myID-1, 2], ICEdata[myID-1, 3], ICEdata[myID-1, 4], ICEdata[myID-1, 5], ICEdata[myID-1, 6], ICEdata[myID-1, 7]
        if  ICEdata[myID-1, 4] == '-':
            prevICE = 'LLAMAICE', ICEdata[myID-2, 2], ICEdata[myID-2, 3], ICEdata[myID-2, 4], ICEdata[myID-2, 5], ICEdata[myID-2, 6], ICEdata[myID-2, 7]
        for i in range(10):
            if (prevICE[1] != '-') & (prevICE[3] != '-'):
                if prevICE[2] in ['-', '']:
                    ib1 = datetime.datetime.strptime(prevICE[1], "%Y-%m-%dT%H:%M" )
                    ib2 = datetime.datetime.strptime(prevICE[3], "%Y-%m-%dT%H:%M" )
                    if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                        axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[0], zorder=0, alpha=alp )
                else:
                    ib1 = datetime.datetime.strptime(prevICE[1], "%Y-%m-%dT%H:%M" )
                    ib2 = datetime.datetime.strptime(prevICE[2], "%Y-%m-%dT%H:%M" )
                    if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                        axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[0], zorder=0, alpha=alp )
                    ib1 = datetime.datetime.strptime(prevICE[2], "%Y-%m-%dT%H:%M" )
                    ib2 = datetime.datetime.strptime(prevICE[3], "%Y-%m-%dT%H:%M" )
                    if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                        axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[1], zorder=0, alpha=alp )    
            elif (prevICE[1] == '-') & (prevICE[2] not in ['-','']):   
                ib1 = datetime.datetime.strptime(prevICE[2], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(prevICE[3], "%Y-%m-%dT%H:%M" )
                axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[1], zorder=0, alpha=alp )
            if (prevICE[3] != '-') & (prevICE[4] != '-'):
                ib1 = datetime.datetime.strptime(prevICE[3], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(prevICE[4], "%Y-%m-%dT%H:%M" )
                if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[2], zorder=0, alpha=alp )
            if (prevICE[4] != '-') & (prevICE[5] != '-'):
                ib1 = datetime.datetime.strptime(prevICE[4], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(prevICE[5], "%Y-%m-%dT%H:%M" )
                if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[3], zorder=0, alpha=alp/2 )
            if (prevICE[5] != '-') & (prevICE[6] != '-'):
                ib1 = datetime.datetime.strptime(prevICE[5], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(prevICE[6], "%Y-%m-%dT%H:%M" )
                if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[4], zorder=0, alpha=alp)
            if (prevICE[5] == '-') & (prevICE[4] != '-') & (prevICE[6] != '-'):
                ib1 = datetime.datetime.strptime(prevICE[4], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(prevICE[6], "%Y-%m-%dT%H:%M" )
                if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[4], zorder=0, alpha=alp)
            # check if sheath overlapping with prev case and shade that region darker
            if (prevICE[-1] != '-') & (ICEbounds[1] != '-'):
                if prevICE[-1] > ICEbounds[1]:
                    ib2 = datetime.datetime.strptime(prevICE[-1], "%Y-%m-%dT%H:%M" )
                    ib1 = datetime.datetime.strptime(ICEbounds[1], "%Y-%m-%dT%H:%M" )
                    if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                        axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[0], zorder=0, alpha=0.1)
                    
        
    # shade in next case
    if (myID != int(totalCMEs)-1) & shadeICE:
        prevICE = 'LLAMAICE', ICEdata[myID+1, 2], ICEdata[myID+1, 3], ICEdata[myID+1, 4], ICEdata[myID+1, 5], ICEdata[myID+1, 6], ICEdata[myID+1, 7]
        if  ICEdata[myID+1, 4] == '-':
            prevICE = 'LLAMAICE', ICEdata[myID+2, 2], ICEdata[myID+2, 3], ICEdata[myID+2, 4], ICEdata[myID+2, 5], ICEdata[myID+2, 6], ICEdata[myID+2, 7]
        for i in range(10):
            if (prevICE[1] != '-') & (prevICE[3] != '-'):
                if prevICE[2] in ['-', '']:
                    ib1 = datetime.datetime.strptime(prevICE[1], "%Y-%m-%dT%H:%M" )
                    ib2 = datetime.datetime.strptime(prevICE[3], "%Y-%m-%dT%H:%M" )
                    if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                        axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[0], zorder=0, alpha=alp )
                else:
                    ib1 = datetime.datetime.strptime(prevICE[1], "%Y-%m-%dT%H:%M" )
                    ib2 = datetime.datetime.strptime(prevICE[2], "%Y-%m-%dT%H:%M" )
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[0], zorder=0, alpha=alp )
                    ib1 = datetime.datetime.strptime(prevICE[2], "%Y-%m-%dT%H:%M" )
                    ib2 = datetime.datetime.strptime(prevICE[3], "%Y-%m-%dT%H:%M" )
                    if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                        axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[1], zorder=0, alpha=alp )    
            elif (prevICE[1] == '-') & (prevICE[2] not in ['-','']):
                ib1 = datetime.datetime.strptime(prevICE[2], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(prevICE[3], "%Y-%m-%dT%H:%M" )
                axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[1], zorder=0, alpha=alp )
            if (prevICE[3] != '-') & (prevICE[4] != '-'):
                ib1 = datetime.datetime.strptime(prevICE[3], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(prevICE[4], "%Y-%m-%dT%H:%M" )
                if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[2], zorder=0, alpha=alp )
            if (prevICE[4] != '-') & (prevICE[5] != '-'):
                ib1 = datetime.datetime.strptime(prevICE[4], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(prevICE[5], "%Y-%m-%dT%H:%M" )
                if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[3], zorder=0, alpha=alp/2 )
            if (prevICE[5] != '-') & (prevICE[6] != '-'):
                ib1 = datetime.datetime.strptime(prevICE[5], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(prevICE[6], "%Y-%m-%dT%H:%M" )
                if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[4], zorder=0, alpha=alp)
            if (prevICE[5] == '-') & (prevICE[4] != '-') & (prevICE[6] != '-'):
                ib1 = datetime.datetime.strptime(prevICE[4], "%Y-%m-%dT%H:%M" )
                ib2 = datetime.datetime.strptime(prevICE[6], "%Y-%m-%dT%H:%M" )
                if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                    axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[4], zorder=0, alpha=alp)
            # check if next sheath overlapping with this case and shade that region darker
            if (ICEbounds[-1] != '-') & (prevICE[1] != '-'):
                if ICEbounds[-1] > prevICE[1]:
                    ib2 = datetime.datetime.strptime(ICEbounds[-1], "%Y-%m-%dT%H:%M" )
                    ib1 = datetime.datetime.strptime(prevICE[1], "%Y-%m-%dT%H:%M" )
                    if (ib1 < endPlotDT) & (ib2 > startPlotDT):
                        axes[i].fill_between([ib1, ib2], bounds[i][0], bounds[i][1], color=cols[4], zorder=0, alpha=0.1)

      
    
# ============================== Add in HSS flags ============================== 
if plotHSS:    
    HSSdata = np.genfromtxt('HSSsort.dat', dtype=str)
    
    yrst = pstart.year
    yrLen = (datetime.datetime(yrst+1,1,1,0,0,0) - datetime.datetime(yrst,1,1,0,0,0)).total_seconds()
    roughstart = yrst + (pstart - datetime.datetime(yrst-1,12,31,0,0,0)).total_seconds() / yrLen

    yren = pend.year
    yrLen = (datetime.datetime(yren+1,1,1,0,0,0) - datetime.datetime(yren,1,1,0,0,0)).total_seconds()
    roughend = yren + (pend - datetime.datetime(yren-1,12,31,0,0,0)).total_seconds() / yrLen

    HSSidxA = np.where((HSSdata[:,1].astype(float) > roughstart) & (HSSdata[:,1].astype(float) < roughend))
    HSSidxB = np.where((HSSdata[:,0].astype(float) > roughstart) & (HSSdata[:,0].astype(float) < roughend))
    HSSidxC = np.where((HSSdata[:,0].astype(float) < roughstart) & (HSSdata[:,1].astype(float) > roughend))
    HSSidx  = np.unique(np.concatenate(HSSidxA+HSSidxB+HSSidxC))
    myHSScats = []
    
    axes[-1].axis('off')
    axes[-1].set_ylim(-0.2, 0.1)
    for idx in HSSidx:
        myCat = HSSdata[idx,2]
        dy = catStyles[myCat][0]
        if myCat not in myHSScats:
            myHSScats.append(myCat)
        if myCat != 'D':
            st_tm = datetime.datetime.strptime(HSSdata[idx,3]+'T'+HSSdata[idx,4], "%Y-%m-%dT%H:%M:%S" )
            en_tm = datetime.datetime.strptime(HSSdata[idx,5]+'T'+HSSdata[idx,6], "%Y-%m-%dT%H:%M:%S" )
        if myCat == 'D':
            st_tm = datetime.datetime.strptime(HSSdata[idx,3]+'T'+HSSdata[idx,4], "%Y-%m-%dT%H:%M:%S" ) 
            en_tm = st_tm + datetime.timedelta(hours=2)           
        axes[-1].fill_between([st_tm, en_tm], [dy-0.1,dy-0.1], [dy,dy], facecolor=catStyles[myCat][1])
           
    
# ============================== Labels and other plotting fun time ============================== 
if addLabels:
    # Catalog Bounds
    nowy = 0.5
    axes[0].text(0.858, nowy,  'Catalog Bounds' , horizontalalignment='left', verticalalignment='center', transform = fig.transFigure, color='k', weight='bold')
    axes[0].annotate('', xy=(0.856, nowy-0.007), xycoords='figure fraction', xytext=(0.9823, nowy-0.007), arrowprops=dict(arrowstyle="-", color='k', lw=1.25))
    nowy -= 0.023
    for i in range(len(printNames)):
        axes[0].text(0.87, nowy, printNames[i], horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color=catCols[printNames[i]], weight='bold')
        nowy -= 0.02
    if ICEdata[myID][1] != '-':
        axes[0].text(0.862, nowy-0.02, '$^*$Nearby Events', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color='gray', weight='bold')
    prevPBend = datetime.datetime.strptime(plotbounds[bidx[0]-1][3], "%Y-%m-%dT%H:%M" )
    nextPBstr = datetime.datetime.strptime(plotbounds[bidx[0]+1][2], "%Y-%m-%dT%H:%M" )
    if (prevPBend > pstart) or (nextPBstr < pend):
        axes[0].text(0.862, nowy-0.02, '$^*$Nearby Events', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color='gray', weight='bold')

    # ICE regions
    nowy = 0.7
    axes[0].text(0.858, nowy, 'ICE Regions', horizontalalignment='left', verticalalignment='center', transform = fig.transFigure, color='k', weight='bold')
    axes[0].annotate('', xy=(0.856, nowy-0.007), xycoords='figure fraction', xytext=(0.9525, nowy-0.007), arrowprops=dict(arrowstyle="-", color='k', lw=1.25))
    nowy -= 0.02
    if ICEcase and shadeICE:
        # Sheath
        if ICEbounds[1] != '-':
            axes[0].text(0.87, nowy, 'Sheath', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color=cols[0], alpha=0.5, weight='bold') 
            nowy -= 0.02
        # Mixed Sheath
        if (ICEbounds[2] != '-') & (ICEbounds[2] != ''):
            axes[0].text(0.87, nowy, 'Mixed', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color=cols[1], alpha=0.5, weight='bold') 
            nowy -= 0.02
        # FR    
        if (ICEbounds[3] != '-'):
            # check if two part or not
            if (ICEbounds[5] == '-'):
                axes[0].text(0.87, nowy, 'FR', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color=cols[2], alpha=0.5, weight='bold') 
                nowy -= 0.02
            else:
                axes[0].text(0.87, nowy, 'FR1', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color=cols[2], alpha=0.5, weight='bold') 
                axes[0].text(0.87, nowy-0.02, 'FR2', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color=cols[3], alpha=0.5/2, weight='bold') 
                nowy -= 0.04
        # mixed at back
        if (ICEbounds[6] != '-'):
            axes[0].text(0.87, nowy, 'Mixed', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color=cols[4], alpha=0.5, weight='bold')
            nowy -= 0.02

        # ICE friends
        friends = ICEdata[myID][1]
        if ';' in friends:
            axes[0].text(0.87, nowy-0.02, '+Prev CME', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color='k', weight='bold')
            axes[0].text(0.87, nowy-0.04, '+Trail CME', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color='k', weight='bold')
        elif friends != '-':
            friend = int(friends)
            if friend < int(CMEchoice):
                axes[0].text(0.87, nowy-0.02, '+Prev CME', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color='k', weight='bold')
            else:
                axes[0].text(0.87, nowy-0.02, '+Trail CME', horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color='k', weight='bold')
            
            
    # add in HSS stuff
    if plotHSS and myHSScats:
        axes[0].text(0.858, 0.97, 'HSS Regions ', horizontalalignment='left', verticalalignment='center', transform = fig.transFigure, color='k', weight='bold')
        nowy = 0.91
        fullnames = {'D':'DONKI', 'V':'VARSITI', 'X':'Xystouris', 'G':'Grandin'}
        axes[0].text(0.858, nowy, 'HSS Catalogs', horizontalalignment='left', verticalalignment='center', transform = fig.transFigure, color='k', weight='bold')
        axes[0].annotate('', xy=(0.856, nowy-0.007), xycoords='figure fraction', xytext=(0.964,nowy-0.007), arrowprops=dict(arrowstyle="-", color='k', lw=1.25))
        nowy -= 0.023
        for cat in myHSScats:
            axes[0].text(0.87, nowy, fullnames[cat], horizontalalignment='left', verticalalignment='center', transform =fig.transFigure, color=catStyles[cat][1], weight='bold')
            nowy -=0.02
        
    
    # Add in obs labels
    axes[0].text(0.858, 0.18, 'Observations', horizontalalignment='left', verticalalignment='center', transform = fig.transFigure, color='k', weight='bold')
    axes[0].annotate('', xy=(0.856, 0.173), xycoords='figure fraction', xytext=(0.964, 0.173), arrowprops=dict(arrowstyle="-", color='k', lw=1.25)) 
    axes[0].text(0.87, 0.157,  'Wind' , horizontalalignment='left', verticalalignment='center', transform = fig.transFigure, color=Windcol, weight='bold')       
    if hasACE:
        axes[0].text(0.87, 0.137,  'ACE' , horizontalalignment='left', verticalalignment='center', transform = fig.transFigure, color=ACEcol, weight='bold')   
            

axes[0].set_ylabel('B [nT]')
axes[1].set_ylabel('Bx [nT]')
axes[2].set_ylabel('By [nT]')
axes[3].set_ylabel('Bz [nT]')
axes[4].set_ylabel('$\\theta$ [$^{\\circ}$]')
axes[5].set_ylabel('$\\phi$ [$^{\\circ}$]')
axes[6].set_ylabel('v [km/s]')
axes[7].set_ylabel('T [K]')
axes[8].set_ylabel('n [cm-3]')
axes[9].set_ylabel('Beta')

axes[9].xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%H:%M"))
#_ = plt.xticks(rotation=45) 
if addLabels:
    plt.subplots_adjust(hspace=0.1,left=0.1,right=0.85,top=0.99,bottom=0.1)
else:
    plt.subplots_adjust(hspace=0.1,left=0.1,right=0.99,top=0.99,bottom=0.1)
axes[9].text(0.035, 0.07, myYr, transform = fig.transFigure)
# add my comments (For now)
if ICEcase:
    axes[9].text(0.05, 0.01, ICEdata[myID][-1], transform = fig.transFigure)    
plt.show()
