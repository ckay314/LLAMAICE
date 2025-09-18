import numpy as np
import datetime
import matplotlib.pyplot as plt
import sys
import matplotlib.dates as mdates

sys.path.append('/Users/kaycd1/LLAMAVERSE/LLAMAICE/FR_fitting')
import lmfit
import prep_data
import rope_fit
import optimise_fit
import output_fit

np.seterr(divide='ignore')

obsDataPath = './ISdata/'

# The other CME catalogs
otherData = np.genfromtxt('AllCases_1.0.dat', dtype=str)
totalCMEs = otherData[-1,0]

# LLAMAICE results which are separate bc more bounds
ICEdata   = np.genfromtxt('LLAMAICE_1.0.csv', delimiter=',', dtype=str)


global ICEcase, failCase
ICEcase = True
failCase = False

global myID, myoIDs, myYr, myDate, myDT, CMEchoice

# Check the number/date given and make sure it's an event we can analyze
def checkCMEchoice(CMEchoice):
    # ============================== ID whether ICE case or generic ==============================
    # =====  check if date or CME number =====
    # assume given a date if has - in it
    global myID, myoIDs, myYr, myDate, myDT, ICEcase
    if '-' in CMEchoice: 
        # should probably make this check columns 3/4 as well for mixed/FR start 
        if CMEchoice in ICEdata[:,2]:
            myID = np.where(ICEdata[:,2] == CMEchoice)[0][0]
            myoIDs = np.where(otherData[:,1] == CMEchoice)[0]
            #print (CMEchoice + ' corresponds to LLAMAICE case '+ICEdata[myID,0])
            # Get the appropriate data files
            myYr = ICEdata[myID,2][0:4]
            myDate = ICEdata[myID,2]
            myDT   = datetime.datetime.strptime(myDate, "%Y-%m-%dT%H:%M" )
        else:
            # This is diff from ICE viewer version bc don't let in non-events here
            sys.exit('Not a LLAMAICE date string so cannot pull bounds and calculate values')            
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
                #print('No LLAMAICE bounds only other catalog(s).')
                dateStr = otherData[myoIDs[0],1]
            myYr = dateStr[0:4]
            myDate = dateStr
            #print ('Results for CME '+dateStr)
            myDT   = datetime.datetime.strptime(dateStr, "%Y-%m-%dT%H:%M" )
        else:
            #print ('Incorrect CME id. LLAMAICE currently contains '+totalCMEs+' CMEs.')
            sys.exit()


def getACEdata(CMEchoice):
    # ===== pull in obs data and format ===== 
    global myYr
    hasACE = True
    if int(myYr) < 1997:
        hasACE = False
    elif int(myYr) == 1997:
        if myDT < datetime.datetime(1997,8,25,0,0,0):
            hasACE = False
    if hasACE:
        aceDataFull = np.genfromtxt(obsDataPath+'ace5/ace_5min_'+myYr+'.dat', dtype=str)
    
    # Select to 2 days before and 2 after, don't need as long as plot range
    startObsDT = myDT - datetime.timedelta(days=3)
    startStr    = startObsDT.strftime("%Y-%m-%dT00:00")  
    endObsDT   = myDT + datetime.timedelta(days=3)
    endStr    = endObsDT.strftime("%Y-%m-%dT23:55")
    # Check if we are going over a year boundary
    if startObsDT.year != myDT.year:
        if int(myYr) >= 1998:
            bonusACE = np.genfromtxt(obsDataPath+'ace5/ace_5min_'+str(int(myYr)-1)+'.dat', dtype=str)
            aceDataFull = np.concatenate((bonusACE, aceDataFull), axis=0)
    elif endObsDT.year != myDT.year:
        if int(myYr) <= 2022:
            bonusACE = np.genfromtxt(obsDataPath+'ace5/ace_5min_'+str(int(myYr)+1)+'.dat', dtype=str)
            aceDataFull = np.concatenate((bonusACE, aceDataFull), axis=0)
    # Extract the data
    try:
        startidx = np.where(aceDataFull[:,0] == startStr)[0][0]
    except:
        startidx = 0
    try:
        endidx  = np.where(aceDataFull[:,0] == endStr)[0][0]
    except:
        endidx = len(aceDataFull[:,0])
    aceData = aceDataFull[startidx:endidx+1, :]
    
    # Format the data
    aceTimes = []
    for i in range(len(aceData)):
        newDT = datetime.datetime.strptime(aceData[i,0], "%Y-%m-%dT%H:%M" )
        aceTimes.append(newDT)
    aceTimes = np.array(aceTimes)
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
    
    return [aceTimes, BxACE, ByACE, BzACE, BtotACE, vsACE, TsACE, nsACE, betaACE, BphiA, BthetaA], [BidxA, vidxA, TidxA, nidxA, betaIdxA]
    
def getWindData(CMEchoice):
    windDataFull = np.genfromtxt(obsDataPath+'wind5/wind_MFISWE_5min_'+myYr+'.dat', dtype=str)
    
    # Select to 2 days before and 2 after, don't need as long as plot range
    startObsDT = myDT - datetime.timedelta(days=3)
    startStr    = startObsDT.strftime("%Y-%m-%dT00:00")  
    endObsDT   = myDT + datetime.timedelta(days=3)
    endStr    = endObsDT.strftime("%Y-%m-%dT23:55")
    # Check if we are going over a year boundary
    if startObsDT.year != myDT.year:
        bonusWind = np.genfromtxt(obsDataPath+'wind5/wind_MFISWE_5min_'+str(int(myYr)-1)+'.dat', dtype=str)
        windDataFull = np.concatenate((windDataFull, bonusWind), axis=0)
    elif endObsDT.year != myDT.year:
        if int(myYr) <= 2022:
            bonusWind = np.genfromtxt(obsDataPath+'wind5/wind_MFISWE_5min_'+str(int(myYr)+1)+'.dat', dtype=str)
            windDataFull = np.concatenate(( windDataFull, bonusWind), axis=0)
    # Extract the data
    try:    
        startidx = np.where(windDataFull[:,0] == startStr)[0][0]
    except:
        startidx = 0 
    try:
        endidx  = np.where(windDataFull[:,0] == endStr)[0][0]
    except:
        endidx = len(windDataFull[:,0])
    windData = windDataFull[startidx:endidx+1, :]
    
    # Format the data
    windTimes = []
    for i in range(len(windData)):
        newDT = datetime.datetime.strptime(windData[i,0], "%Y-%m-%dT%H:%M" )
        windTimes.append(newDT)
    windTimes = np.array(windTimes)
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
    TidxW = np.where((TsWind < 5e6) & (TsWind > 0))[0]
    nidxW = np.where((nsWind < 200) & (nsWind  > 0))[0]
    betaIdxW = []
    for idx in BidxW:
        if idx in TidxW:
            betaIdxW.append(idx)
    # calc plasma beta
    betaWind = (4.16e-5 * TsWind[betaIdxW] + 5.34) * nsWind[betaIdxW]/BtotWind[betaIdxW]**2 # omni version of calc w/Tp, np, dunno where 5.34 comes from
    
    return [windTimes, BxWind, ByWind, BzWind, BtotWind, vsWind, TsWind, nsWind, betaWind, BphiW, BthetaW], [BidxW, vidxW, TidxW, nidxW,betaIdxW]
            
def doThatFit(dataFR, t, theta0, phi0, p0, h, b0, t0):
    Bx, By, Bz = dataFR[0], dataFR[1], dataFR[2] 
    
    fparam = lmfit.Parameters()
    fparam.add('theta0', value=theta0, min=-90, max=90)
    fparam.add('phi0', value=phi0, min=0, max=360)
    fparam.add('p0', value=p0, min=-1, max=1)
    fparam.add('h', value=h, vary=False)
    fparam.add('b0', value=b0, min=0)
    fparam.add('t0', value=t0, min=0)
    
    # Megafit done in one step
    btot_fit, br_fit, bt_fit, bn_fit, vel_fit = rope_fit.fit_initial_guess(fparam, t, dataFR)
    
    # Calculate both chi squared
    chi_dir, chi_mag = optimise_fit.get_chi(Bx, By, Bz, br_fit, bt_fit, bn_fit, b0)
    # Step 1 --- field direction
    solution = lmfit.minimize(optimise_fit.fit_forever, params=fparam, method='leastsq', args=(t,), kws={'data':dataFR, 'keyword':'step1'})

    mparam = lmfit.Parameters()
    mparam.add('theta0', value=solution.params['theta0'].value, vary=False)
    mparam.add('phi0', value=solution.params['phi0'].value, vary=False)
    mparam.add('p0', value=solution.params['p0'].value, vary=False)
    mparam.add('h', value=solution.params['h'].value, vary=False)
    mparam.add('b0', value=solution.params['b0'].value, min=0)
    mparam.add('t0', value=solution.params['t0'].value, vary=False)

    # Step 2 --- field magnitude
    solution = lmfit.minimize(optimise_fit.fit_forever, params=mparam, method='leastsq', args=(t,), kws={'data':dataFR, 'keyword':'step2'})

    newpars = lmfit.create_params(
        theta0 = solution.params['theta0'].value,
        phi0 = solution.params['phi0'].value,
        p0 = solution.params['p0'].value,
        h = solution.params['h'].value,
        b0 = solution.params['b0'].value,
        t0 = solution.params['t0'].value
        )
    
    # Megafit done in one step
    btot_fit, br_fit, bt_fit, bn_fit, vel_fit = rope_fit.fit_initial_guess(newpars, t, dataFR)

    # Calculate both chi squared
    chi_dir, chi_mag = optimise_fit.get_chi(Bx, By, Bz, br_fit, bt_fit, bn_fit, newpars['b0'].value)
    
    return solution, chi_dir, chi_mag

def getRopeDir(t, Bx, By, Bz, v):
    # get avg vectors at front/end of FR to improve init guess
    nT = 12
    bF = [np.mean(Bx[:nT]), np.mean(By[:nT]), np.mean(Bz[:nT])]
    bE = [np.mean(Bx[-nT:]), np.mean(By[-nT:]), np.mean(Bz[-nT:])]
    # mid is harder, could improve using v to find precise mid?
    Tlen = len(Bx)
    mid = int(Tlen / 2)
    bM = [np.mean(Bx[mid-int(nT/2):mid+int(nT/2)+1]), np.mean(By[mid-int(nT/2):mid+int(nT/2)+1]), np.mean(Bz[mid-int(nT/2):mid+int(nT/2)+1])]
        
    # get theta/phi for vector at middle as init guess
    thetaGuess = np.arctan2(bM[2], np.sqrt(bM[0]**2 + bM[1]**2))*180/np.pi + 10
    phiGuess = np.arctan2(bM[1], bM[0])*180/np.pi
    if phiGuess < 0:
        phiGuess += 360.
    
    # just set b0 to magnitude of bM... should be in right ballpark
    b0 = np.sqrt(bM[0]**1 + bM[1]**2 + bM[2]**2)        
    
    # Run through twice with positive and negative h, just keep the best chi sq
    # Also try kicking theta/phi in diff directions
    dThetas = [20,0,-20]
    dPhis   = [30, 0, -30]
    hhhs     = [1, -1]

    nCheck = len(dThetas) * len(dPhis) * 2
    chis   = np.empty(nCheck)
    thetas = np.empty(nCheck)
    phis   = np.empty(nCheck)
    ps     = np.empty(nCheck)
    hs     = np.empty(nCheck)
    bs     = np.empty(nCheck)
    tts    = np.empty(nCheck)
    
    # Set up list with data
    dataFR = [Bx, By, Bz, v]
    
    
    ijk = -1
    for i in range(len(dThetas)):
        for j in range(len(dPhis)):
            for k in range(2):
                ijk += 1
                theta0, phi0, p0, h, b0, t0 = thetaGuess+dThetas[i], phiGuess+dPhis[j], 0.5, hhhs[k], b0, t[-1]
                solution, chi_dir, chi_mag = doThatFit(dataFR, t, theta0, phi0, p0, h, b0, t0)

                chis[ijk]   = chi_dir
                thetas[ijk] = solution.params['theta0'].value
                phis[ijk]   = solution.params['phi0'].value
                hs[ijk]     = solution.params['h'].value
                ps[ijk]     = solution.params['p0'].value
                bs[ijk]     = solution.params['b0'].value
                tts[ijk]    = solution.params['t0'].value
    bestIdx = np.where(chis == np.min(chis))[0]
    bestIdx = bestIdx[0]
    return thetas[bestIdx], phis[bestIdx], ps[bestIdx], hs[bestIdx],  bs[bestIdx], tts[bestIdx], chis[bestIdx]
                
            

def getProps(CMEchoice, sat='both'):    
    
    ACEobs, ACEidx = None, None
    if sat in ['ACE', 'Ace', 'ace', 'BOTH', 'Both', 'both']:
        # order is t, Bx, By, Bz, Btot v, T, n, beta, Bphi, Btheta for obs
        # and B, v, t, n, beta for idx
        if myDT > datetime.datetime(1997,8,25,0,0,0):
            try:
                ACEobs, ACEidx = getACEdata(CMEchoice)
            except:
                print ('Cannot get ACE data')
    
    WindObs, WindIdx = None, None
    if sat in ['WIND', 'Wind', 'wind', 'BOTH', 'Both', 'both']:
        try:
            WindObs, WindIdx = getWindData(CMEchoice)
        except:
            print ('Cannot get Wind data')
            
    # ============================== Read in ICE boundaries ==============================
    toDo = []
    ICEbounds = None
    bonusLines = []
    # Check all lines with this ID number
    for idx in myoIDs:
        myVals = otherData[idx,3], otherData[idx,4], otherData[idx,5], otherData[idx,6] 
        toDo.append(myVals)

    for item in toDo:
        shock = item[1]
        if shock == 'None':
            shock = None
        else:
            shock = datetime.datetime.strptime(shock, "%Y-%m-%dT%H:%M:%S")
        fr1   = datetime.datetime.strptime(item[2], "%Y-%m-%dT%H:%M:%S")
        fr2   = datetime.datetime.strptime(item[3], "%Y-%m-%dT%H:%M:%S")
        
        # Calc ACE properties
        if type(ACEobs) != type(None):
            # shock properties
            outStr = CMEchoice.rjust(3) + ' ' + myDate + item[0].rjust(10) + '  ACE'
            if shock:
                means = np.empty(8)
                for i in range(8):
                    if i < 4:
                        goodIdx = ACEidx[0]
                    else:
                        goodIdx = ACEidx[i-3]
                    myT = ACEobs[0][goodIdx]
                    if i != 7:
                        myVal = ACEobs[i+1][goodIdx] 
                    else:
                        myVal = ACEobs[i+1]
                    idxSh = np.where((myT >= shock) & (myT <= fr1))[0]
                    if len(idxSh) > 0:
                        means[i] = np.mean(myVal[idxSh])
                        outStr += '{:11.2f}'.format(means[i])
                    else:
                        outStr += '       None'            
            else:
                 outStr += '       None       None       None       None       None       None       None       None'

            # FR properties
            means = np.empty(8)
            for i in range(8):
                if i < 4:
                    goodIdx = ACEidx[0]
                else:
                    goodIdx = ACEidx[i-3]
                myT = ACEobs[0][goodIdx]
                if i != 7:
                    myVal = ACEobs[i+1][goodIdx] 
                else:
                    myVal = ACEobs[i+1]
                idxFR = np.where((myT >= fr1) & (myT <= fr2))[0]
                if len(idxFR) > 0:
                    means[i] = np.mean(myVal[idxFR])
                    outStr += '{:11.2f}'.format(means[i])
                else:
                    outStr += '       None'            

            print (outStr)
        
        
        # Calc Wind prpoerties    
        if type(WindObs) != type(None):
            outStr = CMEchoice.rjust(3) + ' ' + myDate +item[0].rjust(10) + ' Wind'
            # shock properties
            if shock:
                means = np.empty(8)
                for i in range(8):
                    if i < 4:
                        goodIdx = WindIdx[0]
                    else:
                        goodIdx = WindIdx[i-3]
                    myT = WindObs[0][goodIdx]
                    if i != 7:
                        myVal = WindObs[i+1][goodIdx] 
                    else:
                        myVal = WindObs[i+1]
                    idxSh = np.where((myT >= shock) & (myT <= fr1))[0]
                    if len(idxSh) > 0:
                        means[i] = np.mean(myVal[idxSh])
                        outStr += '{:11.2f}'.format(means[i])
                    else:
                        outStr += '      None'            
            else:
                 outStr += '       None       None       None       None       None       None       None       None'

            # FR properties
            means = np.empty(8)
            for i in range(8):
                if i < 4:
                    goodIdx = WindIdx[0]
                else:
                    goodIdx = WindIdx[i-3]
                myT = WindObs[0][goodIdx]
                if i != 7:
                    myVal = WindObs[i+1][goodIdx] 
                else:
                    myVal = WindObs[i+1]
                idxFR = np.where((myT >= fr1) & (myT <= fr2))[0]
                if len(idxFR) > 0:
                    means[i] = np.mean(myVal[idxFR])
                    outStr += '{:11.2f}'.format(means[i])
                else:
                    outStr += '       None'  
            print (outStr)   
                    
            # Get Wind orientation       
            '''commonIdx = np.intersect1d(WindIdx[0], WindIdx[1])
            timesDT = WindObs[0][commonIdx]
            
            Bx    = -WindObs[1][commonIdx]
            By    = -WindObs[2][commonIdx]
            Bz    =  WindObs[3][commonIdx]
            v     =  WindObs[5][commonIdx]
            # down select to FR range
            FRidx = np.where((timesDT >= fr1) & (timesDT <= fr2))[0] 
            timesDT = timesDT[FRidx]
            times = np.array([(timesDT[k]-timesDT[0]).total_seconds()/3600. for k in range(len(timesDT))])
            # Need to clean it up once more?
            Bx, By, Bz, v = Bx[FRidx], By[FRidx], Bz[FRidx], v[FRidx]
           
            try:
                theta, phi, p0, h, b0, t0 = getRopeDir(times, Bx, By, Bz, v)
                print (item[0].rjust(12), '{:9.2f}'.format(theta), '{:9.2f}'.format(phi), '{:6.3f}'.format(p0), str(h).rjust(5), '{:9.2f}'.format(b0), '{:9.2f}'.format(t0))
            except:
                print(item[0].rjust(12), '     None', '     None', '  None', ' None', '     None','     None')'''
            
            

def getDirs(CMEchoice, sat='both', fileIn=None):            
    ACEobs, ACEidx = None, None
    if sat in ['ACE', 'Ace', 'ace', 'BOTH', 'Both', 'both']:
        # order is t, Bx, By, Bz, Btot v, T, n, beta, Bphi, Btheta for obs
        # and B, v, t, n, beta for idx
        if myDT > datetime.datetime(1997,8,25,0,0,0):
            try:
                ACEobs, ACEidx = getACEdata(CMEchoice)
            except:
                print ('Cannot get ACE data')
    
    WindObs, WindIdx = None, None
    if sat in ['WIND', 'Wind', 'wind', 'BOTH', 'Both', 'both']:
        try:
            WindObs, WindIdx = getWindData(CMEchoice)
        except:
            print ('Cannot get Wind data')
            
    # ============================== Read in ICE boundaries ==============================
    toDo = []
    ICEbounds = None
    bonusLines = []
    # Check all lines with this ID number
    for idx in myoIDs:
        myVals = otherData[idx,3], otherData[idx,4], otherData[idx,5], otherData[idx,6] 
        toDo.append(myVals)

    for item in toDo:
        fr1   = datetime.datetime.strptime(item[2], "%Y-%m-%dT%H:%M:%S")
        fr2   = datetime.datetime.strptime(item[3], "%Y-%m-%dT%H:%M:%S")
        
        # Calc ACE properties
        if type(ACEobs) != type(None):
            # shock properties
            outStr = CMEchoice.rjust(3) + ' ' + myDate + item[0].rjust(10) + '  ACE'
            commonIdx = np.intersect1d(ACEidx[0], ACEidx[1])
            timesDT = ACEobs[0][commonIdx]
            
            Bx    = -ACEobs[1][commonIdx]
            By    = -ACEobs[2][commonIdx]
            Bz    =  ACEobs[3][commonIdx]
            v     =  ACEobs[5][commonIdx]
            # down select to FR range
            FRidx = np.where((timesDT >= fr1) & (timesDT <= fr2))[0] 
            timesDT = timesDT[FRidx]
            times = np.array([(timesDT[k]-timesDT[0]).total_seconds()/3600. for k in range(len(timesDT))])
            # Need to clean it up once more?
            Bx, By, Bz, v = Bx[FRidx], By[FRidx], Bz[FRidx], v[FRidx]
           
            try:
                theta, phi, p0, h, b0, t0, chi = getRopeDir(times, Bx, By, Bz, v)
                print (outStr, '{:9.2f}'.format(theta), '{:9.2f}'.format(phi), '{:8.3f}'.format(p0), str(h).rjust(5), '{:9.2f}'.format(b0), '{:9.2f}'.format(t0), '{:9.3f}'.format(chi))
                if fileIn:
                    fileIn.write(outStr + '{:9.2f}'.format(theta) + '{:9.2f}'.format(phi)+ '{:8.3f}'.format(p0)+ str(h).rjust(5)+ '{:9.2f}'.format(b0)+ '{:9.2f}'.format(t0)+ '{:9.3f}'.format(chi) + '\n')
            except:
                print(outStr, '     None', '     None', '  None', ' None', '     None','     None','     None')
                
                
            
        # Calc Wind prpoerties    
        if type(WindObs) != type(None):
            outStr = CMEchoice.rjust(3) + ' ' + myDate +item[0].rjust(10) + ' Wind'
            commonIdx = np.intersect1d(WindIdx[0], WindIdx[1])
            timesDT = WindObs[0][commonIdx]
            
            Bx    = -WindObs[1][commonIdx]
            By    = -WindObs[2][commonIdx]
            Bz    =  WindObs[3][commonIdx]
            v     =  WindObs[5][commonIdx]
            # down select to FR range
            FRidx = np.where((timesDT >= fr1) & (timesDT <= fr2))[0] 
            timesDT = timesDT[FRidx]
            times = np.array([(timesDT[k]-timesDT[0]).total_seconds()/3600. for k in range(len(timesDT))])
            # Need to clean it up once more?
            Bx, By, Bz, v = Bx[FRidx], By[FRidx], Bz[FRidx], v[FRidx]
           
            try:
                theta, phi, p0, h, b0, t0, chi = getRopeDir(times, Bx, By, Bz, v)
                print (outStr, '{:9.2f}'.format(theta), '{:9.2f}'.format(phi), '{:8.3f}'.format(p0), str(h).rjust(5), '{:9.2f}'.format(b0), '{:9.2f}'.format(t0), '{:9.3f}'.format(chi))
                if fileIn:
                    fileIn.write(outStr + '{:9.2f}'.format(theta) + '{:9.2f}'.format(phi)+ '{:8.3f}'.format(p0)+ str(h).rjust(5)+ '{:9.2f}'.format(b0)+ '{:9.2f}'.format(t0)+ '{:9.3f}'.format(chi) + '\n')
            except:
                print(outStr, '     None', '     None', '  None', ' None', '     None','     None','     None')
                
                
        
        


def runIt():
    global CMEchoice, myDate
    CMEchoice = sys.argv[1]
    if CMEchoice in ['All', 'all', 'ALL']:
        stdout_bak = sys.stdout
        f1 = open('tempOut.txt', 'w')
        #with open('tempOut.txt', 'w') as sys.stdout:
        for CMEchoice in ICEdata[:,0]:
            if int(CMEchoice) > -1: # convient to only loop through part if catches for some reason
                checkCMEchoice(CMEchoice)
                #getProps(CMEchoice)
                getDirs(CMEchoice, fileIn=f1)
        f1.close()
    else:
        f1 = open('tempOut.txt', 'w')
        checkCMEchoice(CMEchoice)
        #getProps(CMEchoice, sat='both')
        getDirs(CMEchoice, fileIn=f1)
        f1.close
        
if __name__ == '__main__':
    runIt()