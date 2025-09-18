import numpy as np


global data, unNums
data = np.genfromtxt('meanData.dat', dtype=str)
# columns are num date user sat sheath[bx by bz B v T n beta] (4-11) fr[bx by bz B v T n beta] (12-19)
unNums = np.unique(data[:,0])
nTot   = len(unNums)

unDate = np.unique(data[:,1])

def getUncs():
    for ii in range(nTot):
        i = ii + 1
        myACE = np.where((data[:,0].astype(int) == i) & (data[:,3] == 'ACE'))[0]
        myWind = np.where((data[:,0].astype(int) == i) & (data[:,3] == 'Wind'))[0]
        outPrintA = ''
        outPrintW = ''
        for jj in range(16):
            j = jj + 4
            
            if i >= 41: # hardcoded for when ACE data starts
                subValsA = data[myACE,j]
                goodsA = np.where(subValsA != 'None')[0]
                subValsA = subValsA[goodsA] 
                if len(subValsA) == 0:
                   if j  in [9,17]:
                       myOutA = '       None' +  '      None'
                   elif j in [8, 16]:
                       myOutA = '     None' +  '   None'
                   else:
                       myOutA = '   None' +  '   None'
                else:
                    if j  in [9,17]:
                        myOutA = '{:11.2f}'.format(np.mean(subValsA.astype(float))) +  '{:10.2f}'.format(np.std(subValsA.astype(float)))
                    elif j in [8, 16]:
                        myOutA = '{:9.2f}'.format(np.mean(subValsA.astype(float))) +  '{:7.2f}'.format(np.std(subValsA.astype(float)))
                    else:
                        myOutA = '{:7.2f}'.format(np.mean(subValsA.astype(float))) +  '{:7.2f}'.format(np.std(subValsA.astype(float)))
                
                outPrintA += myOutA


            subValsW = data[myWind,j]
            goodsW = np.where(subValsW != 'None')[0]
            subValsW = subValsW[goodsW] 
            if len(subValsW) == 0:
               if j  in [9,17]:
                   myOutW = '       None' +  '      None'
               elif j in [8, 16]:
                   myOutW = '     None' +  '   None'
               else:
                   myOutW = '   None' +  '   None'
            else:
                if j  in [9,17]:
                    myOutW = '{:11.2f}'.format(np.mean(subValsW.astype(float))) +  '{:10.2f}'.format(np.std(subValsW.astype(float)))
                elif j in [8, 16]:
                    myOutW = '{:9.2f}'.format(np.mean(subValsW.astype(float))) +  '{:7.2f}'.format(np.std(subValsW.astype(float)))
                else:
                    myOutW = '{:7.2f}'.format(np.mean(subValsW.astype(float))) +  '{:7.2f}'.format(np.std(subValsW.astype(float)))
                
            outPrintW += myOutW
        if i >= 41:
            print (str(i).rjust(3), unDate[ii], '  ACE', outPrintA)    
        print (str(i).rjust(3), unDate[ii], ' Wind', outPrintW)
        #print (sd)

def getStats():
    data1 = np.genfromtxt('singleMeans.dat', dtype = str)
    names = ['Bx (nT)', 'By (nT)', 'Bz (nT)', 'B (nT)', 'v (km/s)', 'T (K)', 'n (cm^-3)', 'Beta']
    for i in range(16):
        idA = 2*i +4
        idB = idA - 1
        # Wind Data
        idx = np.where((data1[:,2] == 'Wind') & (data1[:,idA] != 'None') & (data1[:,idA] != '0.00'))
        subValsW = data1[idx,idB].astype(float)
        subUncsW = data1[idx,idA].astype(float)
        percsW = subUncsW / np.abs(subValsW)
        percsW = percsW[np.isfinite(percsW)]
        
        # ACE Data
        # Wind Data
        idx = np.where((data1[:,2] == 'ACE') & (data1[:,idA] != 'None') & (data1[:,idA] != '0.00'))
        subValsA = data1[idx,idB].astype(float)
        subUncsA = data1[idx,idA].astype(float)
        percsA = subUncsA / np.abs(subValsA)
        percsA = percsA[np.isfinite(percsA)]
        
        if i == 0:
            print( 'Sheath Uncertainties:       ACE                  Wind')
        elif i ==8:
            print ('')
            print( 'FR Uncertainties:           ACE                  Wind')
        print ((names[i%8]+':').rjust(10)+'{:12.2f}'.format(np.mean(subUncsA)), '{:8.2f}'.format(np.mean(percsA)*100), '{:12.2f}'.format(np.mean(subUncsW)), '{:8.2f}'.format(np.mean(percsW)*100))
        
        
            
        
        
        
                

#getUncs()   
getStats() 