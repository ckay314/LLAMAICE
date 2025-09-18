import matplotlib.pyplot as plt
import numpy as np
import datetime 
from matplotlib.ticker import FixedLocator, MultipleLocator
import matplotlib.gridspec as gridspec
from scipy import stats
from sklearn.metrics import r2_score
from scipy.interpolate import interp1d
from scipy.interpolate import make_interp_spline



# make label text size bigger
plt.rcParams.update({'font.size':14})



catCols = {'ATSB':'Cyan', 'CDAW':'Tomato', 'DONKI':'Maroon', 'DREAMS':'Magenta', 'GMU1':'#F0E68C', 'GMU2':'Gold', 'H4C':'#9932CC', 'Jian':'Maroon', 'Lepping':'#F39C12', 'RC':'Blue', 'Wind':'Red', 'LLAMAICEA':'k', 'LLAMAICEB':'RebeccaPurple', 'LLAMAICEC':'Turquoise','LLAMAICED':'PowderBlue'}
catOrds = { 'CDAW':6,  'DREAMS':1,  'Jian':2, 'Lepping':3, 'RC':5, 'Wind':4, 'LLAMAICED':0}

shcol = '#882255'
frcol = '#88CCEE'

def timelinePlot():
    catOrds = { 'CDAW':6,  'DREAMS':5,  'Jian':4, 'Lepping':3, 'RC':2, 'Wind':1, 'LLAMAICEA':0}
    
    fig = plt.figure(figsize=(10,10))
    gs = fig.add_gridspec(12, 1, hspace=0.0, top=0.99, right=0.99, left=0.05, bottom=0.07)
    fs = []
    for i in range(12):
        if i == 0:
            fs.append(fig.add_subplot(gs[i]))
            fs[0].set_yticks([])
            fs[0].tick_params(labelbottom=False)
        else:
            fs.append(fig.add_subplot(gs[i], sharey=fs[0], sharex=fs[0]))
            if i != 11:
                fs[i].tick_params(labelbottom=False)
        fs[i].set_ylabel(1995+i)
    fs[-1].set_xlabel('Day of Year')
    fs[0].set_ylim([-1, 7])
    fs[0].set_xlim([0,367])
    
    data = np.genfromtxt('AllCases_1.0.dat', dtype=str)
    # print out total counts
    for name in catOrds:
        print (name, len(np.where(data[:,3] == name)[0]))
    for i in range(len(data[:,0])):
        myCol, myOrd = None, None
        if data[i, 3] in catOrds:
            myCol = catCols[data[i, 3]]
            myOrd = catOrds[data[i,3]]
        '''elif data[i, 3] == 'LLAMAICEB':
            myCol = catCols['LLAMAICED']
            myOrd = catOrds['LLAMAICED']
        elif data[i, 3] == 'LLAMAICEA':
            myCol = catCols['LLAMAICED']
            myOrd = catOrds['LLAMAICED']'''
        
        
        if data[i,4] != 'None':              
            front = datetime.datetime.strptime(data[i, 4], "%Y-%m-%dT%H:%M:%S" )
        else:
            front = datetime.datetime.strptime(data[i, 5], "%Y-%m-%dT%H:%M:%S" )
        back = datetime.datetime.strptime(data[i, 6], "%Y-%m-%dT%H:%M:%S" )
                    
        myYr = front.year
        doyFront = (front - datetime.datetime(myYr, 1,1)).total_seconds()/3600/24 + 1
        doyBack = (back - datetime.datetime(myYr, 1,1)).total_seconds()/3600/24 + 1
        if doyBack > 366:
            yrMax = datetime.datetime(myYr, 12,31,23, 59, 59)
            if back > yrMax:
                bonus = (back-yrMax).total_seconds()/3600/24 + 1
                fs[int(myYr-1995)+1].plot([1,bonus], [myOrd, myOrd], c=myCol, lw=5)
                
        fs[int(myYr-1995)].plot([doyFront,doyBack], [myOrd, myOrd], c=myCol, lw=5)    
    plt.savefig('paperFigs/timeline.png')
    
def compAW():
    fig, ax = plt.subplots(4,2, figsize=(5.5,10))
    
    
    #gs = gridspec.GridSpec(nrows=4, ncols=2)
    #gs = fig.add_gridspec(4, 2, hspace=0.2, wspace=0.2, top=0.99, right=0.99, left=0.1, bottom=0.1)
    fs = []
    for k in range(8):
        i,j = int(k/2), k%2
        #fs.append(fig.add_subplot(gs[i,j]))
        ax[i,j].set_aspect('equal')
        
    
    
    data = np.genfromtxt('singleMeans.dat', dtype=str)
    nCases = int(data[-1,0])
    shIdx = [3,5,7,9,11,13,15,17]
    frIdx = [19,21,23,25,27,29,31,33]
    aceSh = [[] for i in range(8)]
    windSh = [[] for i in range(8)]
    aceFR = [[] for i in range(8)]
    windFR = [[] for i in range(8)]
    for i in range(nCases):
        j = i+1
        idx = np.where(data[:,0] == str(j))[0]
        if len(idx) == 2:
            for k in range(8):
                ks = shIdx[k]
                if (data[idx[0], ks] != 'None') & (data[idx[1], ks] != 'None'):
                    aceSh[k].append(float(data[idx[1], ks]))
                    windSh[k].append(float(data[idx[0], ks]))    
                kfr = frIdx[k]
                if (data[idx[0], kfr] != 'None') & (data[idx[1], kfr] != 'None'):
                    aceFR[k].append(float(data[idx[1], kfr]))
                    windFR[k].append(float(data[idx[0], kfr]))
    for i in range(8):
        aceFR[i] = np.array(aceFR[i])                             
        windFR[i] = np.array(windFR[i])                             
        aceSh[i] = np.array(aceSh[i])                             
        windSh[i] = np.array(windSh[i])                             
    for i in [5,7]:     
        windSh[i] = np.log10(windSh[i])                             
        aceSh[i] = np.log10(aceSh[i])                             
        windFR[i] = np.log10(windFR[i])                             
        aceFR[i] = np.log10(aceFR[i])                             
    myLabs = ['Bx (nT)', 'By (nT)', 'Bz (nT)', 'B (nT)', 'v (km/s)', 'log T (K)', 'n (cm$^{-3}$)', 'log Beta']    
    setRg = {0:[-25,25], 1:[-25,25], 2:[-25,25],3:[0,50], 6:[0,70]}            
    for k in range(8):
        i,j = int(k/2), k%2
        ax[i,j].scatter(aceSh[k], windSh[k], c=shcol, s=10, marker='x')
        ax[i,j].scatter(aceFR[k], windFR[k], c=frcol, s=10, marker='x')
        print(myLabs[k])
        if k != 7:
            #print(stats.ttest_ind(aceSh[k], windSh[k]))
            #print(stats.ttest_ind(aceFR[k], windFR[k]))
            _, pnormS = stats.mannwhitneyu(aceSh[k], windSh[k], use_continuity=False,    method="auto")
            _, pnormF = stats.mannwhitneyu(aceFR[k], windFR[k], use_continuity=False,    method="auto")
            
            aS, aF = aceSh[k], aceFR[k]
            wS, wF = windSh[k], windFR[k]
            r, p = stats.pearsonr(aS, wS)
            print('{:.3f}'.format(pnormS), '{:.2f}'.format(np.mean(np.abs(aS - wS))),'{:.3f}'.format(r), '{:.2e}'.format(p),   '{:.3f}'.format(r2_score(wS,aS)), '{:.3f}'.format(r2_score(aS, wS)))
            r, p = stats.pearsonr(aF, wF)
            print('{:.3f}'.format(pnormF),'{:.2f}'.format(np.mean(np.abs(aF - wF))),  '{:.3f}'.format(r), '{:.2e}'.format(p),  '{:.3f}'.format(r2_score(wF,aF)), '{:.3f}'.format(r2_score(aF, wF)))
            #print(np.mean(np.abs(aceSh[k] - windSh[k])), stats.pearsonr(aceSh[k], windSh[k]), pnormS, r2_score(aceSh[k], windSh[k]), r2_score(windSh[k], aceSh[k]))
            #print(np.mean(np.abs(aceFR[k] - windFR[k])), stats.pearsonr(aceFR[k], windFR[k]), pnormF, r2_score(aceFR[k], windFR[k]), r2_score(windFR[k], aceFR[k]))
            
            
        if k ==7:
            aS, aF = np.power(10,aceSh[k]), np.power(10,aceFR[k])
            wS, wF =  np.power(10,windSh[k]), np.power(10,windFR[k])
            _, pnormS = stats.mannwhitneyu(aS, wS, use_continuity=False,    method="auto")
            _, pnormF = stats.mannwhitneyu(aF, wF, use_continuity=False,    method="auto")
            #print(stats.ttest_ind(aS,wS), pnorm)
            #print(stats.ttest_ind(aF, wS))
            r, p = stats.pearsonr(aS, wS)
            print('{:.3f}'.format(pnormS), '{:.2f}'.format(np.mean(np.abs(aS - wS))), '{:.3f}'.format(r), '{:.2e}'.format(p),  '{:.3f}'.format(r2_score(wS,aS)), '{:.3f}'.format(r2_score(aS, wS)))
            r, p = stats.pearsonr(aF, wF)
            print('{:.3f}'.format(pnormF), '{:.2f}'.format(np.mean(np.abs(aF - wF))), '{:.3f}'.format(r), '{:.2e}'.format(p),  '{:.3f}'.format(r2_score(wF,aF)), '{:.3f}'.format(r2_score(aF, wF)))
            #print(np.mean(np.abs(aF - wF)), stats.pearsonr(aF, wF), pnormF, r2_score(aF, wF), r2_score(wF, aF))
        
        if k in setRg:
            newrng = setRg[k]
        else:
            xr, yr = ax[i,j].get_xlim(), ax[i,j].get_ylim()
            newrng = [np.min([xr[0], yr[0]]), np.max([xr[1], yr[1]])]
        ax[i,j].plot(newrng, newrng, 'k--', zorder=0)
        ax[i,j].set_xlim(newrng)
        ax[i,j].set_ylim(newrng)
        ax[i,j].set_xlabel('ACE '+myLabs[k])
        ax[i,j].set_ylabel('Wind '+myLabs[k])

    fig.subplots_adjust( hspace=0.5, wspace=0.1, top=0.99, right=0.99, left=0.1, bottom=0.1)            
    #plt.show()
    plt.savefig('paperFigs/AWscatter.png')

def ratioPlot():
    data = np.genfromtxt('meanData.dat', dtype=str)
    nCases = int(data[-1,0])
    
    names = np.unique(data[:,2])
    outSh = {}
    outFR = {}
    for name in names:
        if name != 'LLAMAICEA':
            outSh[name] = [[] for i in range(8)]
            outFR[name] = [[] for i in range(8)]
            
    shIdx = [4,5,6,7,8,9,10,11]
    frIdx = [12,13,14,15,16,17,18,19]    
    for i in range(nCases):
        j = i+1
        idx = np.where((data[:,0] == str(j)) & (data[:,3] == 'Wind'))[0]
        myNames = data[idx,2]
        if 'LLAMAICEA' in myNames:
            LAid = idx[np.where(myNames == 'LLAMAICEA')[0]][0]
            idx = np.delete(idx, np.where(myNames == 'LLAMAICEA')[0])
            for k in idx:
                for l in range(8):
                    if (data[LAid, shIdx[l]] != 'None') & (data[k, shIdx[l]] != 'None'):
                        thisRat = float(data[k, shIdx[l]]) / float(data[LAid, shIdx[l]])
                        outSh[data[k,2]][l].append(thisRat)
                    if (data[LAid, frIdx[l]] != 'None') & (data[k, frIdx[l]] != 'None'):
                        thisRat = float(data[k, frIdx[l]]) / float(data[LAid, frIdx[l]])
                        outFR[data[k,2]][l].append(thisRat)
    



    name2ordSh = {'CDAW':0, 'DREAMS':1, 'Jian':2, 'RC':3, 'Wind':4, 'LLAMAICEC':5}
    name2ordFR = {'CDAW':0, 'DREAMS':1, 'Jian':2, 'Lepping':3, 'RC':4, 'Wind':5, 'LLAMAICEB':6, 'LLAMAICEC':7, 'LLAMAICED':8}                  
    shOuts = np.zeros([5,6])
    frOuts = np.zeros([5,9])  
    for key in outSh.keys():
        if key in name2ordSh:
            for j in range(5):
                i = j+3
                print (key, 'Sh', len(outSh[key][i]), np.mean(outSh[key][i]))
                shOuts[j,name2ordSh[key]] = (np.mean(outSh[key][i]) - 1) * 100
        if key in name2ordFR:
            for j in range(5):
                i = j+3
                print (key, 'FR', len(outFR[key][i]), np.mean(outFR[key][i]))
                frOuts[j,name2ordFR[key]] = (np.mean(outFR[key][i]) - 1) * 100
                
    fig = plt.figure(figsize=(8.5, 3.5))
    gs = fig.add_gridspec(1, 10)
    f1 = fig.add_subplot(gs[0,0:4]) 
    f2 = fig.add_subplot(gs[0,4:9], sharey=f1) 
    f3 = fig.add_subplot(gs[0,9])
    
    vval = 50
    f1.pcolor(shOuts[::-1,:], vmin=-vval, vmax=vval, cmap='seismic')
    im = f2.pcolor(frOuts[::-1,:], vmin=-vval, vmax=vval, cmap='seismic')
    cb = fig.colorbar(im, cax=f3)
    
    # Making pretty
    cb.set_label('Percent Variation', rotation=270, labelpad=20)
    #f2.set_yticks([])
    f1.set_yticks([0.5 + i for i in range(5)],['Beta', 'n', 'T', 'v', 'B'])
    plt.setp(f2.get_yticklabels(), visible=False)
    xlabSh = ['CDAW', 'DREAMS', 'Jian', 'RC', 'Wind', 'LLAMA$_C$']
    xlabFR = ['CDAW', 'DREAMS', 'Jian', 'Lepping', 'RC', 'Wind', 'LLAMA$_B$', 'LLAMA$_C$', 'LLAMA$_D$']

    f1.set_xticks([0.5 + i for i in range(6)], xlabSh, fontsize=10, rotation =270)
    f2.set_xticks([0.5 + i for i in range(9)], xlabFR, fontsize=10, rotation =270)
    f1.set_title('Sheath')
    f2.set_title('Flux Rope')
    fig.subplots_adjust(  top=0.85, right=0.9, left=0.1, bottom=0.25)            
    
    plt.savefig('paperFigs/ratios.png')

def catHistos():
    data = np.genfromtxt('meanData.dat', dtype=str)
    nCases = int(data[-1,0])
    satName = 'Wind'
    
    names = np.unique(data[:,2])
    outSh = {}
    outFR = {}
    for name in names:
        outSh[name] = [[] for i in range(8)]
        outFR[name] = [[] for i in range(8)]
            
    shIdx = [4,5,6,7,8,9,10,11]
    frIdx = [12,13,14,15,16,17,18,19]    
    for i in range(nCases):
        j = i+1
        idx = np.where((data[:,0] == str(j)) & (data[:,3] == satName))[0]
        myNames = data[idx,2]
        for k in idx:
                for l in range(8):
                    if (data[k, shIdx[l]] != 'None'):
                        thisRat = float(data[k, shIdx[l]]) 
                        outSh[data[k,2]][l].append(thisRat)
                    if  (data[k, frIdx[l]] != 'None'):
                        thisRat = float(data[k, frIdx[l]])
                        outFR[data[k,2]][l].append(thisRat)
    



    name2ordSh = {'CDAW':0, 'DREAMS':1, 'Jian':2, 'RC':3, 'Wind':4, 'LLAMAICEA':5, 'LLAMAICEC':6}
    name2ordFR = {'CDAW':0, 'DREAMS':1, 'Jian':2, 'Lepping':3, 'RC':4, 'Wind':5, 'LLAMAICEA':6, 'LLAMAICEB':7, 'LLAMAICEC':8, 'LLAMAICED':9}     
    
    fig, axes = plt.subplots(4,2, figsize=(7, 9))
    fs = [axes[0,0], axes[0,1], axes[1,0], axes[1,1], axes[2,0], axes[2,1], axes[3,0], axes[3,1]]


    Bval = 25
    binRng = [[-Bval,Bval], [-Bval,Bval], [-Bval,Bval], [0,35], [300,800], [4,6.5], [0,40], [-1.5,1.5]]
    myLabs = ['Bx (nT)', 'By (nT)', 'Bz (nT)', 'B (nT)', 'v (km/s)', 'log T (K)', 'n (cm$^{-3}$)', 'log Beta']
    myLabs2 = ['Bx', 'By', 'Bz', 'B', 'v', 'log T', 'n', 'log Beta']

    # FR plot
    handles = []
    for i in range(8):
        ps = []
        ordKeys = []
        for key in ['CDAW', 'DREAMS', 'Jian', 'Lepping', 'RC', 'Wind', 'LLAMAICEA', 'LLAMAICEB', 'LLAMAICEC', 'LLAMAICED']:
            myVals = outFR[key][i]
            if i in [5,7]:
                myVals = np.log10(myVals)
            ae, loce, scalee = stats.skewnorm.fit(myVals)
            x = np.linspace(binRng[i][0], binRng[i][-1], 100)
            p = stats.skewnorm.pdf(x, ae, loce, scalee)
            ps.append(p)
            ordKeys.append(key)
        
        pmax = np.max(ps)
        for ii in range(len(ps)):
            p = ps[ii] / pmax
            x = np.linspace(binRng[i][0], binRng[i][-1], 100)
            if i == 1:
                line, = fs[i].plot(x, p, '--', c=catCols[ordKeys[ii]], linewidth=3, label=ordKeys[ii])
                handles.append(line)
            else:
                fs[i].plot(x, p, '--', c=catCols[ordKeys[ii]], linewidth=3)
    for i in range(8):
        fs[i].set_xlim(binRng[i][0], binRng[i][-1])
        fs[i].set_xlabel(myLabs[i])
        if i % 2 == 0:
            fs[i].set_ylabel('ScPDF')
    
    labels = [h.get_label() for h in handles]
    fig.legend(handles, labels, loc='upper center', ncol=4, fontsize=11)      
    fig.subplots_adjust(hspace=0.5, wspace=0.25,  top=0.9, right=0.98, left=0.15, bottom=0.07)    
    plt.savefig('paperFigs/catFRpdf_'+satName+'.png')     
    #plt.show()
    #print (outSh.keys())    
    
    # Sheath plot
    fig, axes = plt.subplots(4,2, figsize=(7, 9))
    fs = [axes[0,0], axes[0,1], axes[1,0], axes[1,1], axes[2,0], axes[2,1], axes[3,0], axes[3,1]]
    handles = []   
    for i in range(8):
        ps = []
        ordKeys = []
        for key in ['CDAW', 'DREAMS', 'Jian', 'RC', 'Wind', 'LLAMAICEA', 'LLAMAICEC']:
            myVals = np.array(outSh[key][i])
            if i in [5,7]:
                myVals = np.log10(myVals)
            myVals = myVals[np.isfinite(myVals)]
            ae, loce, scalee = stats.skewnorm.fit(myVals)
            x = np.linspace(binRng[i][0], binRng[i][-1], 100)
            p = stats.skewnorm.pdf(x, ae, loce, scalee)    
            ps.append(p)
            ordKeys.append(key)        
        
        pmax = np.max(ps)    
        for ii in range(len(ps)):
            p = ps[ii] / pmax    
            if i == 1:
                line, = fs[i].plot(x, p, '--', c=catCols[ordKeys[ii]], linewidth=3, label=ordKeys[ii])
                handles.append(line)
            else:
                fs[i].plot(x, p, '--', c=catCols[ordKeys[ii]], linewidth=3)
    for i in range(8):
        fs[i].set_xlim(binRng[i][0], binRng[i][-1])
        fs[i].set_xlabel(myLabs[i])
        if i % 2 == 0:
            fs[i].set_ylabel('ScPDF')
    
    labels = [h.get_label() for h in handles]
    fig.legend(handles, labels, loc='upper center', ncol=4, fontsize=11)      
    fig.subplots_adjust(hspace=0.5, wspace=0.25,  top=0.9, right=0.98, left=0.15, bottom=0.07)    
    #plt.show()
    plt.savefig('paperFigs/catShpdf_'+satName+'.png')
    
    allStats = np.zeros([8,10,10])*np.nan
    
    # Print out all the stats
    print('Sheath MWU stats: ')
    for i in range(8):
        print(myLabs[i])
        for key1 in ['CDAW', 'DREAMS', 'Jian', 'LLAMAICEA', 'LLAMAICEC', 'RC', 'Wind']:
            for key2 in ['CDAW', 'DREAMS', 'Jian', 'LLAMAICEA', 'LLAMAICEC', 'RC', 'Wind']:
                if key2 != key1:
                    _, pnormS = stats.mannwhitneyu(outSh[key1][i], outSh[key2][i], use_continuity=False,    method="auto")
                    print(key1, key2, np.log10(pnormS))
                    idx1, idx2 = name2ordFR[key1], name2ordFR[key2] 
                    if idx1 > idx2:
                        allStats[i, idx1, idx2 ] = np.log10(pnormS)
        print ('')
    
    print ('')
    print ('')
    print ('')
    
    print('FR MWU stats: ')
    for i in range(8):
        print(myLabs[i])
        for key1 in outFR:
            for key2 in outFR:
                if key2 != key1:
                    _, pnormS = stats.mannwhitneyu(outFR[key1][i], outFR[key2][i], use_continuity=False,    method="auto")
                    print(key1, key2, np.log10(pnormS))
                    idx1, idx2 = name2ordFR[key1], name2ordFR[key2] 
                    print(key1, key2, np.log10(pnormS), idx1+1, idx2+1)
                    if idx1 < idx2:
                        allStats[i, idx1, idx2 ] = np.log10(pnormS)
        print ('')
    
    fig = plt.figure(figsize=(5.5, 9))
    gs = fig.add_gridspec(5, 2, height_ratios=[1,1, 1,1,0.1])
    fs = []
    for i in range(4):
        fs.append(fig.add_subplot(gs[i,0]))
        fs.append(fig.add_subplot(gs[i,1]))
    cax = fig.add_subplot(gs[4,:])
    #fs = [axes[0,0], axes[0,1], axes[1,0], axes[1,1], axes[2,0], axes[2,1], axes[3,0], axes[3,1]]
    for i in range(8):
        if i == 0:
            im = fs[i].pcolormesh(allStats[i,:9,:], vmin=-5, vmax=0, cmap='plasma_r',edgecolors='k',lw=0.7)
        else:
            fs[i].pcolormesh(allStats[i,:9,:], vmin=-5, vmax=0, cmap='plasma_r',edgecolors='k', lw=0.7)
        fs[i].set_aspect('equal')
        fs[i].set_xticks([0.5,2.5,4.5,6.5,8.5],labels=[1,3,5,7,9])
        fs[i].set_yticks([1.5,3.5,5.5,7.5],labels=[2,4,6,8])
        fs[i].set_title(myLabs2[i],fontsize=14)
        
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    #cax.yaxis.set_ticks_position('left')
    cax.set_title('MWU p',fontsize=14)  
    cb.set_ticks([0,-1,-2,-3,-4,-5])
    cb.set_ticklabels(['1', '10$^{-1}$', '10$^{-2}$','10$^{-3}$','10$^{-4}$','10$^{-5}$']) 
    fig.subplots_adjust(hspace=0.5, wspace=0.1,  top=0.95, right=0.95, left=0.05, bottom=0.07)  
    plt.savefig('paperFigs/MWU'+satName+'.png')
    
                
def paramHistos():
    satTag = 'ACE'
    
    data = np.genfromtxt('singleMeans.dat', dtype=str)
    shIdx = [3,5,7,9,11,13,15,17] 
    frIdx = [19,21,23,25,27,29,31,33] 
    
    # reusing uncertainty code and not bothering to change param names
    # but these are params/means not uncertainties
    ASunc = [[] for i in range(8)]
    AFunc = [[] for i in range(8)]
    WSunc = [[] for i in range(8)]
    WFunc = [[] for i in range(8)]
    for i in range(len(data[:,0])):
        myL = data[i,:]
        if myL[2] == 'ACE':
            for j in range(8):
                if myL[shIdx[j]+1] not in ['None']:
                    ASunc[j].append(float(myL[shIdx[j]])) 
                if myL[frIdx[j]+1] not in ['None']:
                    AFunc[j].append(float(myL[frIdx[j]])) 
        elif myL[2] == 'Wind':
            for j in range(8):
                if myL[shIdx[j]+1] not in ['None', '0.00']:
                    WSunc[j].append(float(myL[shIdx[j]])) 
                if myL[frIdx[j]+1] not in ['None', '0.00']:
                    WFunc[j].append(float(myL[frIdx[j]])) 
        
    fig, ax = plt.subplots(8,2, figsize=(6.5,10))
    fsh = ax[:,0]
    ffr = ax[:,1]
    
    for i in range(8):
        ASunc[i] = np.array(ASunc[i])
        AFunc[i] = np.array(AFunc[i])
        WSunc[i] = np.array(WSunc[i])
        WFunc[i] = np.array(WFunc[i])
    
    ASunc[5] = np.log10(ASunc[5])
    AFunc[5] = np.log10(AFunc[5])
    WSunc[5] = np.log10(WSunc[5])
    WFunc[5] = np.log10(WFunc[5]) 
    
    ASunc[7] = np.log10(ASunc[7])
    AFunc[7] = np.log10(AFunc[7])
    WSunc[7] = np.log10(WSunc[7])
    WFunc[7] = np.log10(WFunc[7]) 
    
    Bval = 25
    binRng = [[-Bval,Bval], [-Bval,Bval], [-Bval,Bval], [0,35], [300,800], [4,6.5], [0,40], [-1.5,1.5]]
    myLabs = ['Bx (nT)', 'By (nT)', 'Bz (nT)', 'B (nT)', 'v (km/s)', 'log T (K)', 'n (cm$^{-3}$)', 'log Beta']    
    if satTag == 'ACE':
        myS = ASunc
        myF = AFunc
        maxC = 160
    elif satTag == 'Wind':
        myS = WSunc
        myF = WFunc
        maxC = 140
    
    
    for i in range(8):
        bins = np.linspace(binRng[i][0],binRng[i][1],10 )
        '''myS[i][np.where(myS[i] > bins[-1])] = 0.5 * (bins[-2] + bins[-1])
        myF[i][np.where(myF[i] > bins[-1])] = 0.5 * (bins[-2] + bins[-1])
        
        myS[i][np.where(myS[i] < bins[0])] = 0.5 * (bins[0] + bins[1])
        myF[i][np.where(myF[i] < bins[0])] = 0.5 * (bins[0] + bins[1])'''
        
        ae, loce, scalee = stats.skewnorm.fit(WSunc[i])
        x = np.linspace(binRng[i][0], binRng[i][-1], 100)
        pW = stats.skewnorm.pdf(x, ae, loce, scalee)
        ae, loce, scalee = stats.skewnorm.fit(ASunc[i])
        x = np.linspace(binRng[i][0], binRng[i][-1], 100)
        pA = stats.skewnorm.pdf(x, ae, loce, scalee)
        maxp = np.max([pW, pA])
        fsh[i].plot(x,pW/maxp, 'k--')
        fsh[i].plot(x,pA/maxp, '--', c='gray')
    
        ae, loce, scalee = stats.skewnorm.fit(WFunc[i])
        x = np.linspace(binRng[i][0], binRng[i][-1], 100)
        pW = stats.skewnorm.pdf(x, ae, loce, scalee)
        ae, loce, scalee = stats.skewnorm.fit(AFunc[i])
        x = np.linspace(binRng[i][0], binRng[i][-1], 100)
        pA = stats.skewnorm.pdf(x, ae, loce, scalee)
        maxp = np.max([pW, pA])
        ffr[i].plot(x,pW/maxp, 'k--')
        ffr[i].plot(x,pA/maxp, '--', c='gray')

        fsh[i].set_xlim(binRng[i][0], binRng[i][-1])
        ffr[i].set_xlim(binRng[i][0], binRng[i][-1])
        fsh[i].set_ylim(0, 1.1)
        ffr[i].set_ylim(0, 1.1)
        
        '''fsh[i].hist(myS[i], fc=shcol, ec='k', bins=bins)
        ffr[i].hist(myF[i], fc=frcol, ec='k', bins=bins)
        fsh[i].set_ylim(0,maxC)
        ffr[i].set_ylim(0,maxC)'''
        if i != 7:
            fsh[i].text(0.98, 0.7, myLabs[i], fontsize=12, transform=fsh[i].transAxes, horizontalalignment='right')
        else:
            fsh[i].text(0.05, 0.7, myLabs[i], fontsize=12, transform=fsh[i].transAxes, horizontalalignment='left')
        ffr[i].text(0.98, 0.7, myLabs[i], fontsize=12, transform=ffr[i].transAxes, horizontalalignment='right')
        fsh[i].set_ylabel('ScPDF')
        #plt.setp(ffr[i].get_yticklabels(), visible=False)
    ffr[0].set_title('Flux Rope')
    fsh[0].set_title('Sheath')
    
    
                
    fig.subplots_adjust(hspace=0.5, wspace=0.2,  top=0.95, right=0.97, left=0.13, bottom=0.05)            
    plt.savefig('paperFigs/paramHistosAW.png')
    #plt.show()

def uncHistos():
    satTag = 'ACE'
    
    data = np.genfromtxt('singleMeans.dat', dtype=str)
    shIdx = [3,5,7,9,11,13,15,17] # these are mean vals not unc
    frIdx = [19,21,23,25,27,29,31,33] # unc is just +1
    
    ASunc = [[] for i in range(8)]
    AFunc = [[] for i in range(8)]
    WSunc = [[] for i in range(8)]
    WFunc = [[] for i in range(8)]
    for i in range(len(data[:,0])):
        myL = data[i,:]
        if myL[2] == 'ACE':
            for j in range(8):
                if myL[shIdx[j]+1] not in ['None', '0.00']:
                    ASunc[j].append(float(myL[shIdx[j]+1])) 
                if myL[frIdx[j]+1] not in ['None', '0.00']:
                    AFunc[j].append(float(myL[frIdx[j]+1])) 
        elif myL[2] == 'Wind':
            for j in range(8):
                if myL[shIdx[j]+1] not in ['None', '0.00']:
                    WSunc[j].append(float(myL[shIdx[j]+1])) 
                if myL[frIdx[j]+1] not in ['None', '0.00']:
                    WFunc[j].append(float(myL[frIdx[j]+1])) 
        
    fig, ax = plt.subplots(8,2, figsize=(5,10))
    fsh = ax[:,0]
    ffr = ax[:,1]
    
    for i in range(8):
        ASunc[i] = np.array(ASunc[i])
        AFunc[i] = np.array(AFunc[i])
        WSunc[i] = np.array(WSunc[i])
        WFunc[i] = np.array(WFunc[i])
    
    ASunc[5] = ASunc[5] / 1e4
    AFunc[5] = AFunc[5] / 1e4
    WSunc[5] = WSunc[5] / 1e4
    WFunc[5] = WFunc[5] / 1e4
    
    binRng = [[0,3], [0,4], [0,4], [0,3], [0,30], [0,5], [0,4], [0,0.75]]
    myLabs = ['Bx (nT)', 'By (nT)', 'Bz (nT)', 'B (nT)', 'v (km/s)', 'T (10$^4$ K)', 'n (cm$^{-3}$)', 'Beta']    
    if satTag == 'ACE':
        myS = ASunc
        myF = AFunc
        maxC = 120
    elif satTag == 'Wind':
        myS = WSunc
        myF = WFunc
        maxC = 140
    
    for i in range(8):
        
        bins = np.linspace(binRng[i][0],binRng[i][1],10 )
        myS[i][np.where(myS[i] > bins[-1])] = 0.5 * (bins[-2] + bins[-1])
        myF[i][np.where(myF[i] > bins[-1])] = 0.5 * (bins[-2] + bins[-1])
        fsh[i].hist(myS[i], bins=bins, fc=shcol, ec='k')
        ffr[i].hist(myF[i], bins=bins, fc=frcol, ec='k')
        fsh[i].set_ylim(0,maxC)
        ffr[i].set_ylim(0,maxC)
        fsh[i].text(0.95, 0.7, myLabs[i], fontsize=12, transform=fsh[i].transAxes, horizontalalignment='right')
        ffr[i].text(0.95, 0.7, myLabs[i], fontsize=12, transform=ffr[i].transAxes, horizontalalignment='right')
        fsh[i].set_ylabel('Counts')
        plt.setp(ffr[i].get_yticklabels(), visible=False)
    ffr[0].set_title('Flux Rope')
    fsh[0].set_title('Shock')
    
    
                
    fig.subplots_adjust(hspace=0.5, wspace=0.1,  top=0.95, right=0.99, left=0.15, bottom=0.05)            
    plt.savefig('paperFigs/uncHistos_'+satTag+'.png')
    #plt.show()
    
    
def fitStats():
    data = np.genfromtxt('FRfits.dat', dtype=str)
    # Check for questionable fits
    counterA = 0
    matchA   = 0
    matchA1  = 0
    counterW = 0
    matchW   = 0
    matchW1  = 0
    
    matchIncA = 0
    matchIdxA = np.zeros(3)
    matchIncA1 = 0
    matchIdxA1 = np.zeros(3)
    matchIncW = 0
    matchIdxW = np.zeros(3)
    matchIncW1 = 0
    matchIdxW1 = np.zeros(3)
    
    AWmatches = 0
    AWtots    = 0
    Aincs = np.zeros(3)
    Wincs = np.zeros(3)
    incMatches = 0
    matchCodes = np.zeros(3)
    
    for i in range(int(data[-1,0])):
        j = i+1
        subIdx = np.where(data[:,0].astype(int) == j )[0]
        subACE = np.where((data[:,0].astype(int) == j) & (data[:,3] == 'ACE') )[0]
        subWind = np.where((data[:,0].astype(int) == j) & (data[:,3] == 'Wind') )[0]
        
        # Found 8 cases where only one fit and LLAMAICE just uses the same as that one
        # 17 cases where gets stuck on 0/180/360 and 2 case where all same but rando vals
        
        
        if len(subACE) > 1:
            if np.std(data[subACE,4].astype(float)) != 0:
                #print ('ACE', j, len(np.unique(data[subACE, 7])))
                counterA += 1
                # Handedness grouping
                if len(np.unique(data[subACE, 7])) == 1:
                    matchA += 1
                else:
                    nums = [len(np.where(data[subACE, 7].astype(float) == 1)[0]), len(np.where(data[subACE, 7].astype(float) == -1)[0])]
                    if (np.min(nums) == 1) & (len(subACE) >2):
                        matchA1 += 1
                # Inclination
                incs = data[subACE,4].astype(float)
                angs = data[subACE,5].astype(float)
                nTot = len(subACE)
                # Pulled these inclination defs from Erika 2018
                nHigh = len(np.where(np.abs(incs) >= 55)[0])
                nLow  = len(np.where(np.abs(incs) <= 35)[0])
                nMid  =  nTot-nHigh-nLow
                ns    = np.array([nLow, nMid, nHigh])
                n0s = len(np.where(ns == 0)[0])
                if n0s == 2:
                    matchIncA += 1
                    idx = np.where(ns !=0)[0]
                    matchIdxA[idx] += 1
                elif (nTot > 2) & (n0s ==1):
                    if np.sort(ns)[1] == 1:
                        idx1 = np.where(ns == np.max(ns))[0]
                        matchIncA1 += 1
                        matchIdxA1[idx1] += 1
                
                        
                        
        
        if len(subWind) > 1:
             if np.std(data[subWind,4].astype(float)) != 0:
                 #print ('Win', j, len(np.unique(data[subACE, 7])))
                 counterW += 1
                 # Handedness grouping
                 if len(np.unique(data[subWind, 7])) == 1:
                     matchW += 1
                 else:
                     nums = [len(np.where(data[subWind, 7].astype(float) == 1)[0]), len(np.where(data[subWind, 7].astype(float) == -1)[0])]
                     if (np.min(nums) == 1) & (len(subWind) > 2):
                         matchW1 += 1
                 # Inclination
                 incs = data[subWind,4].astype(float)
                 angs = data[subWind,5].astype(float)
                 nTot = len(subWind)
                 nHigh = len(np.where(np.abs(incs) >= 55)[0])
                 nLow  = len(np.where(np.abs(incs) <= 35)[0])
                 nMid  =  nTot-nHigh-nLow
                 ns    = np.array([nLow, nMid, nHigh])
                 n0s = len(np.where(ns == 0)[0])
                 if n0s == 2:
                    matchIncW += 1
                    idx = np.where(ns !=0)[0]
                    matchIdxW[idx] += 1
                 elif (nTot > 2) & (n0s ==1):
                    if np.sort(ns)[1] == 1:
                        idx1 = np.where(ns == np.max(ns))[0]
                        matchIncW1 += 1
                        matchIdxW1[idx1] += 1
                        
        # has both ace and wind
        if (len(subACE) >= 1) & (len(subWind) >= 1):
            for idx in subACE:
                idx2 = idx + 1
                if idx2 in subWind:
                #if data[idx, 2] in data[subWind,2]:
                    AWtots += 1
                    #idx2 = subWind[np.where(data[subWind,2] == data[idx,2])[0]]
                    if (data[idx,7] == data[idx2,7]):
                        AWmatches += 1
                    inc1 = float(data[idx,4])
                    inc2 = float(data[idx2,4])
                    # ACE
                    if inc1 <=35:
                        Aincs[0] += 1
                        Acode = 'L'
                    elif inc1 >=55:
                        Aincs[2] += 1
                        Acode = 'H'
                    else:
                        Aincs[1] += 1
                        Acode = 'M'
                    # Wind
                    if inc2 <=35:
                        Wincs[0] += 1
                        Wcode = 'L'
                    elif inc2 >=55:
                        Wincs[2] += 1
                        Wcode = 'H'
                    else:
                        Wincs[1] += 1
                        Wcode = 'M'
                    if Acode == Wcode:
                        incMatches += 1
                        if Acode == 'L':
                            matchCodes[0] += 1
                        elif Acode == 'M':
                            matchCodes[1] += 1
                        else:
                            matchCodes[2] += 1
                    
    
    print ('ACE/Wind results')   
    print ('Handedness match: ', AWmatches, ' out of ', AWtots )
    print ('                  ', AWmatches/ AWtots )
    print ('Inclination match: ', incMatches, ' out of ', AWtots)
    print ('                   ', incMatches / AWtots)
    print ('Matches : ', matchCodes, '[low, mid, high]')
    print (' ACE inc: ', Aincs, '[low, mid, high]')
    print ('          ', Aincs/AWtots)
    print (' Wind inc:', Wincs, '[low, mid, high]')
    print ('          ', Wincs/AWtots)
    
    print ('')              
                     
         
    print ('ACE polarity results')
    print (counterA, ' cases with unique reconstructions')
    print (matchA, ' cases with consistent handedness')
    print ('    ', matchA / counterA)
    print (matchA1, ' cases with one differing')
    print ('    ', matchA1 / counterA)
    print ('')
    print ('Wind polarity results')
    print (counterW, ' cases with unique reconstructions')
    print (matchW, ' cases with consistent handedness')
    print ('    ', matchW / counterW)
    print (matchW1, ' cases with one differing')
    print ('    ', matchW1 / counterW)
                   
    print ('')
    print ('')
    print ('')
    print ('ACE Inclination results')
    print (matchIncA, matchIdxA, 'consistent cases. [low, mid, high]')
    print ('   ', matchIncA/counterA)
    print (matchIncA1, matchIdxA1, 'one off cases')
    print ('   ', matchIncA1/counterA)
    print ('')
    print ('Wind Inclination results')
    print (matchIncW, matchIdxW, 'consistent cases. [low, mid, high]')
    print ('   ', matchIncW/counterW)
    print (matchIncW1, matchIdxW1, 'one off cases')
    print ('   ', matchIncW1/counterW)
    
def compFitAW():
    data = np.genfromtxt('FRfits.dat', dtype=str)
    
    allIncDiffs = []
    allIncAs = []
    allIncWs = []
    IncCols = []

    allLonDiffs = []
    allLonAs = []
    allLonWs = []
    LonCols = []
    
    counter = 0
    for i in range(int(data[-1,0])+1):
        j = i+1
        idx = np.where(data[:,0].astype(int) == j)
        subIdx = np.where(data[:,0].astype(int) == j )[0]
        subACE = np.where((data[:,0].astype(int) == j) & (data[:,3] == 'ACE') )[0]
        subWind = np.where((data[:,0].astype(int) == j) & (data[:,3] == 'Wind') )[0]
        
        #names = np.unique(data[subIdx,2])
        if (len(subACE) >= 1) & (len(subWind) >= 1):
            for aceIdx in subACE:
                windIdx = aceIdx + 1
                if windIdx in subWind:
 
                    counter += 1

                    name = data[aceIdx,2]
                    aceInc = float(data[aceIdx,4])
                    windInc = float(data[windIdx,4])
                    diff = np.abs(aceInc - windInc)
                    
                    if diff > 90:
                        diff = 180 - diff
                    
                    if True:
                    #if np.abs(aceInc) + np.abs(windInc) != 0:
                        allIncDiffs.append(diff)
                        allIncAs.append(aceInc)
                        allIncWs.append(windInc)
                        if name in catCols:
                            IncCols.append(catCols[name])
                        else:
                            IncCols.append(catCols['LLAMAICED'])
                    
                    aceLon = float(data[aceIdx,5])
                    windLon = float(data[windIdx,5])
                    diff = np.abs(aceLon - windLon)
                    if diff > 180:
                        diff = 360 - diff
                    
                    if True:
                    #if (aceLon != windLon):# and (aceLon not in [0, 180, 360]):
                        allLonDiffs.append(diff)
                        allLonAs.append(aceLon)
                        if (windLon - aceLon) > 180:
                            windLon -= 360
                        elif (windLon - aceLon) < -180:
                            windLon += 360
                        allLonWs.append(windLon)
                        if name in catCols:
                            LonCols.append(catCols[name])
                        else:
                            LonCols.append(catCols['LLAMAICED'])
                        
    allIncDiffs = np.array(allIncDiffs)
    allLonDiffs = np.array(allLonDiffs)
                    
    fig = plt.figure(constrained_layout=True, figsize=(12, 8.5))
    gs = fig.add_gridspec(3, 2)
    f1 = fig.add_subplot(gs[0:2,0]) 
    f2 = fig.add_subplot(gs[0:2,1]) 
    f3a = fig.add_subplot(gs[2:,0]) 
    f3b = fig.add_subplot(gs[2:,1], sharey=f3a) 
    
    
    f1.scatter(allIncAs, allIncWs, c=IncCols)                
    f2.scatter(allLonAs, allLonWs, c=LonCols) 
    
    bmax1 = 20
    bins1 = np.linspace(0,bmax1, 20)
    allIncDiffs[np.where(allIncDiffs[i] > bins1[-1])] = 0.5 * (bins1[-2] + bins1[-1])
    f3a.hist(allIncDiffs, bins=bins1, ec='k')               

    bmax2 = 40
    bins2 = np.linspace(0,bmax2, 20)
    allLonDiffs[np.where(allLonDiffs[i] > bins2[-1])] = 0.5 * (bins2[-2] + bins2[-1])
    f3b.hist(allLonDiffs, bins=bins2, ec='k')               
 
    f1.plot([-90,90], [-90,90], 'k--', zorder=0) 
    f2.plot([-40,400], [-40,400], 'k--', zorder=0) 
 
    f1.set_xlim([-90,90])
    f1.set_ylim([-90,90])
    f2.set_xlim([-40,400])
    f2.set_ylim([-40,400])
    f1.set_aspect('equal')
    f2.set_aspect('equal')
    
    f3a.set_xticks([0,4,8,12,16,20])
    f3a.set_xlim([0,bmax1])
    f3b.set_xlim([0,bmax2])
    
    f1.set_xlabel('ACE $\\theta$ ($^{\\circ}$)')
    f1.set_ylabel('Wind $\\theta$ ($^{\\circ}$)')
    f2.set_xlabel('ACE $\\phi$ ($^{\\circ}$)')
    f2.set_ylabel('Wind $\\phi$ ($^{\\circ}$)')
    
    f3a.set_xlabel('|$\\Delta \\theta$| ($^{\\circ}$)')
    f3b.set_xlabel('|$\\Delta \\phi$| ($^{\\circ}$)')
    f3a.set_ylabel('Counts')
    f3b.set_ylabel('Counts')

    plt.savefig('paperFigs/fitAW.png')     
    print (counter)
    print ('Inc 5,10 ', len(np.where(allIncDiffs < 5)[0]), len(np.where(allIncDiffs < 10)[0]), len(allIncDiffs))  
    print (' ', len(np.where(allIncDiffs < 5)[0])/len(allIncDiffs), len(np.where(allIncDiffs < 10)[0])/len(allIncDiffs) )       
    print ('Lon <5,10 ', len(np.where(allLonDiffs < 5)[0]), len(np.where(allLonDiffs < 10)[0])) 
    print (' ', len(np.where(allLonDiffs < 5)[0])/len(allLonDiffs), len(np.where(allLonDiffs < 10)[0])/len(allLonDiffs))    
    print ("")   
    print (np.mean(allIncDiffs), np.median(allIncDiffs), np.std(allIncDiffs), len(allIncDiffs))
    print (np.mean(allLonDiffs), np.median(allLonDiffs), np.std(allLonDiffs), len(allLonDiffs))

def compLLAMA():    
    data = np.genfromtxt('FRfits.dat', dtype=str)
    satTag = 'ACE'
    
    thetas = [[] for i in range(4)]
    phis = [[] for i in range(4)]

    for i in range(int(data[-1,0])):
        j = i+1
        idx = np.where(data[:,0].astype(int) == j)
        subIdx = np.where((data[:,0].astype(int) == j) & (data[:,3]==satTag) )[0]
        names = np.unique(data[subIdx,2])
        isLlama = ['LLAMA' in name for name in names]
        llamaNames = names[isLlama]
        # Find the multi-Llama cases
        if (len(llamaNames) > 1) & ('LLAMAICEA' in llamaNames):
            thisCase = data[subIdx,:]
            # LA
            thisIdx = np.where((thisCase[:,2] == llamaNames[0]) & (thisCase[:,3]==satTag))[0]
            thetas[0].append(thisCase[thisIdx[0],4])
            phis[0].append(thisCase[thisIdx[0],5])
            # LB
            if 'LLAMAICEB' in llamaNames:
                thisIdx = np.where((thisCase[:,2] == 'LLAMAICEB') & (thisCase[:,3]==satTag))[0]
                thetas[1].append(thisCase[thisIdx[0],4])
                phis[1].append(thisCase[thisIdx[0],5])
            else:
                thetas[1].append(None)
                phis[1].append(None)
            # LC
            if 'LLAMAICEC' in llamaNames:
                thisIdx = np.where((thisCase[:,2] == 'LLAMAICEC') & (thisCase[:,3]==satTag))[0]
                thetas[2].append(thisCase[thisIdx[0],4])
                phis[2].append(thisCase[thisIdx[0],5])
            else:
                thetas[2].append(None)
                phis[2].append(None)
            
            # LD
            if 'LLAMAICED' in llamaNames:
                thisIdx = np.where((thisCase[:,2] == 'LLAMAICED') & (thisCase[:,3]==satTag))[0]
                thetas[3].append(thisCase[thisIdx[0],4])
                phis[3].append(thisCase[thisIdx[0],5])
            else:
                thetas[3].append(None)
                phis[3].append(None)
    
    
    for i in range(4):
        thetas[i] = np.array(thetas[i]).astype(float)
        phis[i] = np.array(phis[i]).astype(float)
    
    diffsTheta = [[] for i in range(3)]
    diffsPhi = [[] for i in range(3)]
    # Clean up angles        
    for i in range(len(phis[0])):
        for jj in range(3):
            j = jj+1
            if phis[j][i] - phis[0][i]  > 180:
                phis[j][i] -=360
            elif phis[j][i] - phis[0][i]  < -180:
                phis[j][i] += 360
            #if j == 3:
            #    print (phis[j][i], phis[0][i])
            
            if np.isfinite(thetas[j][i]):
                dTheta = np.abs(thetas[j][i] - thetas[0][i])
                if dTheta > 90:
                    dTheta = 180 - dTheta
                diffsTheta[jj].append(dTheta)
            
            if np.isfinite(phis[j][i]):
                dPhi = np.abs(phis[j][i] - phis[0][i])
                if dPhi > 180:
                    dPhi = 360 - dPhi
                diffsPhi[jj].append(dPhi)
                
            if phis[j][i] < -80:
                phis[j][i] += 360
            if phis[j][i] > 400:
                phis[j][i] -=360
                
    fig, ax = plt.subplots(3,2, figsize=(7,10))
    for i in range(3):
        ax[i,0].scatter(thetas[0], thetas[i+1])
        ax[i,1].scatter(phis[0], phis[i+1])
        ax[i,0].plot([-90, 90], [-90,90], 'k--', zorder=0)
        ax[i,1].plot([-60,400], [-60,400], 'k--', zorder=0)
        
        ax[i,0].set_aspect('equal')
        ax[i,1].set_aspect('equal')
        ax[i,0].set_xlim([-90,90])
        ax[i,0].set_ylim([-90,90])
        ax[i,1].set_xlim([-60,400])
        ax[i,1].set_ylim([-60,400])
        
        ax[i,0].set_xlabel('$\\theta_A$ ($^{\\circ}$)')
        ax[i,1].set_xlabel('$\\phi_A$ ($^{\\circ}$)')
        
    ax[0,0].set_ylabel('$\\theta_B$ ($^{\\circ}$)')
    ax[0,1].set_ylabel('$\\phi_B$ ($^{\\circ}$)')
    ax[1,0].set_ylabel('$\\theta_C$ ($^{\\circ}$)')
    ax[1,1].set_ylabel('$\\phi_C$ ($^{\\circ}$)')
    ax[2,0].set_ylabel('$\\theta_D$ ($^{\\circ}$)')
    ax[2,1].set_ylabel('$\\phi_D$ ($^{\\circ}$)')

    fig.subplots_adjust(hspace=0.3, wspace=0.1,  top=0.95, right=0.99, left=0.1, bottom=0.1)            
    
    plt.savefig('paperFigs/LLAMAfits_'+satTag+'.png')
    for i in range(3):
        print (i, np.mean(diffsTheta[i]), np.mean(diffsPhi[i]), np.median(diffsTheta[i]), np.median(diffsPhi[i]), len(diffsTheta[i]))

def compRecons():
    data = np.genfromtxt('FRfits.dat', dtype=str)
    
    thetaW  = []
    thetaUW = []
    thetaA  = []
    thetaUA = []
    phiW    = []
    phiUW   = []
    phiA    = []
    phiUA   = []
    for i in range(int(data[-1,0])):
        j = i + 1
        idxW = np.where((data[:,0].astype(int) == j) & (data[:,3] == 'Wind'))[0]
        idxA = np.where((data[:,0].astype(int) == j) & (data[:,3] == 'ACE'))[0]
        if len(np.unique(data[idxW,4].astype(float))) > 1:
            # Inclination (the easy one)
            myIncs = data[idxW,4].astype(float)
            thetaW.append(np.mean(myIncs))
            thetaUW.append(np.std(myIncs))
            
            myLons = data[idxW,5].astype(float)
            if (np.max(myLons)-np.min(myLons)) > 180.:
                myLons[np.where(myLons > 180)] -= 360.
            phiW.append(np.mean(myLons))
            phiUW.append(np.std(myLons))
            
        if len(np.unique(data[idxA,4].astype(float))) > 1:
            # Inclination (the easy one)
            myIncs = data[idxA,4].astype(float)
            thetaA.append(np.mean(myIncs))
            thetaUA.append(np.std(myIncs))
            
            myLons = data[idxA,5].astype(float)
            if (np.max(myLons)-np.min(myLons)) > 180.:
                myLons[np.where(myLons > 180)] -= 360.
            phiA.append(np.mean(myLons))
            phiUA.append(np.std(myLons))
            
    thetaW  = np.array(thetaW)
    thetaUW = np.array(thetaUW)
    thetaA  = np.array(thetaA)
    thetaUA = np.array(thetaUA)
    phiW    = np.array(phiW)
    phiUW   = np.array(phiUW)
    phiA    = np.array(phiA)
    phiUA   = np.array(phiUA)      
      
    print ('Theta Wind:', np.mean(thetaUW), np.median(thetaUW))
    print ('Theta ACE: ', np.mean(thetaUA), np.median(thetaUA))
    print ('Phi Wind:  ', np.mean(phiUW), np.median(phiUW))
    print ('PHI ACE:   ', np.mean(phiUA), np.median(phiUA))
    
    fig, ax = plt.subplots(2,2, figsize=(7,7))
    bmaxes = [75, 150]
    bins1 = np.linspace(0,bmaxes[0], 15)
    bins2 = np.linspace(0,bmaxes[1], 15)
    ax[0,0].hist(thetaUW, bins=bins1, ec='k')
    ax[1,0].hist(thetaUA, bins=bins1, ec='k')
    ax[0,1].hist(phiUW, bins=bins2, ec='k')
    ax[1,1].hist(phiUA, bins=bins2, ec='k')
    
    for i in [0,1]:
        for j in [0,1]:
            ax[i,j].set_ylabel('Counts')
            ax[i,j].set_xlim([0,bmaxes[j]])
    ax[0,0].set_xlabel('$\\Delta\\theta_{Wind}$')
    ax[1,0].set_xlabel('$\\Delta\\theta_{ACE}$')
    ax[0,1].set_xlabel('$\\Delta\\phi_{Wind}$')
    ax[1,1].set_xlabel('$\\Delta\\phi_{ACE}$')
            
    
    fig.subplots_adjust(hspace=0.3, wspace=0.3,  top=0.95, right=0.97, left=0.1, bottom=0.1)            
    
    plt.savefig('paperFigs/fitUnc.png')

def compLLAMArecons():
    data = np.genfromtxt('FRfits.dat', dtype=str)
    satTag = 'Wind'
    
    thetas = [[] for i in range(4)]
    phis = [[] for i in range(4)]

    for i in range(int(data[-1,0])):
        j = i+1
        idx = np.where(data[:,0].astype(int) == j)
        subIdx = np.where((data[:,0].astype(int) == j) & (data[:,3]==satTag) )[0]
        names = np.unique(data[subIdx,2])
        isLlama = ['LLAMA' in name for name in names]
        llamaNames = names[isLlama]
        # Find the multi-Llama cases
        if (len(llamaNames) > 1) & ('LLAMAICEA' in llamaNames):
            thisCase = data[subIdx,:]
            # LA
            thisIdx = np.where((thisCase[:,2] == llamaNames[0]) & (thisCase[:,3]==satTag))[0]
            thetas[0].append(thisCase[thisIdx[0],4])
            phis[0].append(thisCase[thisIdx[0],5])
            # LB
            if 'LLAMAICEB' in llamaNames:
                thisIdx = np.where((thisCase[:,2] == 'LLAMAICEB') & (thisCase[:,3]==satTag))[0]
                thetas[1].append(thisCase[thisIdx[0],4])
                phis[1].append(thisCase[thisIdx[0],5])
            else:
                thetas[1].append(None)
                phis[1].append(None)
            # LC
            if 'LLAMAICEC' in llamaNames:
                thisIdx = np.where((thisCase[:,2] == 'LLAMAICEC') & (thisCase[:,3]==satTag))[0]
                thetas[2].append(thisCase[thisIdx[0],4])
                phis[2].append(thisCase[thisIdx[0],5])
            else:
                thetas[2].append(None)
                phis[2].append(None)
            
            # LD
            if 'LLAMAICED' in llamaNames:
                thisIdx = np.where((thisCase[:,2] == 'LLAMAICED') & (thisCase[:,3]==satTag))[0]
                thetas[3].append(thisCase[thisIdx[0],4])
                phis[3].append(thisCase[thisIdx[0],5])
            else:
                thetas[3].append(None)
                phis[3].append(None)
    
    
    for i in range(4):
        thetas[i] = np.array(thetas[i]).astype(float)
        phis[i] = np.array(phis[i]).astype(float)
    
    diffsTheta = [[] for i in range(3)]
    diffsPhi = [[] for i in range(3)]
    # Clean up angles        
    for i in range(len(phis[0])):
        for jj in range(3):
            j = jj+1
            if phis[j][i] - phis[0][i]  > 180:
                phis[j][i] -=360
            elif phis[j][i] - phis[0][i]  < -180:
                phis[j][i] += 360
            #if j == 3:
            #    print (phis[j][i], phis[0][i])
            
            if np.isfinite(thetas[j][i]):
                dTheta = np.abs(thetas[j][i] - thetas[0][i])
                if dTheta > 90:
                    dTheta = 180 - dTheta
                diffsTheta[jj].append(dTheta)
            
            if np.isfinite(phis[j][i]):
                dPhi = np.abs(phis[j][i] - phis[0][i])
                if dPhi > 180:
                    dPhi = 360 - dPhi
                diffsPhi[jj].append(dPhi)
                
            if phis[j][i] < -80:
                phis[j][i] += 360
            if phis[j][i] > 400:
                phis[j][i] -=360
                
    fig, ax = plt.subplots(3,2, figsize=(7,7))
    bmaxes = [75, 150]
    bins1 = np.linspace(0,bmaxes[0], 15)
    bins2 = np.linspace(0,bmaxes[1], 15)
    for i in range(3):
        ax[i,0].hist(diffsTheta[i], bins=bins1, ec='k',density=True)
        ax[i,1].hist(diffsPhi[i], bins=bins2, ec='k', density=True)
    
    ymaxes = [0.1, 0.05]
    tits = ['B - A', 'C - A', 'D - A']
    for i in [0,1,2]:
        for j in [0,1]:
            ax[i,j].set_xlim([0,bmaxes[j]])
            ax[i,j].set_ylim([0,ymaxes[j]])
            ax[i,j].text(0.95, 0.8, tits[i], transform=ax[i,j].transAxes, horizontalalignment='right')
        ax[i,0].set_ylabel('Prob. Dens.')
    ax[2,0].set_xlabel('$\\Delta\\theta$')
    ax[2,1].set_xlabel('$\\Delta\\phi$')
    
    
            
    
    fig.subplots_adjust(hspace=0.3, wspace=0.3,  top=0.95, right=0.97, left=0.15, bottom=0.1)            
    #plt.show()
    plt.savefig('paperFigs/fitUncLLAMA'+satTag+'.png')
    

def fitPDFs():
    data = np.genfromtxt('FRfits.dat', dtype=str)
    
    
    # Compare all recons for ACE/Wind
    allIncDiffs = []
    allIncAs = []
    allIncWs = []

    allLonDiffs = []
    allLonAs = []
    allLonWs = []

    for i in range(int(data[-1,0])+1):
        j = i+1
        idx = np.where(data[:,0].astype(int) == j)
        subIdx = np.where(data[:,0].astype(int) == j )[0]
        subACE = np.where((data[:,0].astype(int) == j) & (data[:,3] == 'ACE') )[0]
        subWind = np.where((data[:,0].astype(int) == j) & (data[:,3] == 'Wind') )[0]
        
        if (len(subACE) >= 1) & (len(subWind) >= 1):
            for aceIdx in subACE:
                windIdx = aceIdx + 1
                if windIdx in subWind:
                    name = data[aceIdx,2]
                    aceInc = float(data[aceIdx,4])
                    windInc = float(data[windIdx,4])
                    diff = np.abs(aceInc - windInc)
                    
                
                    if diff > 90:
                        diff = 180 - diff
                
                    if np.abs(aceInc) + np.abs(windInc) != 0:
                        allIncDiffs.append(diff)
                        allIncAs.append(aceInc)
                        allIncWs.append(windInc)
                                           
                    aceLon = float(data[aceIdx,5])
                    windLon = float(data[windIdx,5])
                    diff = np.abs(aceLon - windLon)
                    if diff > 180:
                        diff = 360 - diff
                    
                    if (aceLon != windLon):# and (aceLon not in [0, 180, 360]):
                        allLonDiffs.append(diff)
                        allLonAs.append(aceLon)
                        if (windLon - aceLon) > 180:
                            windLon -= 360
                        elif (windLon - aceLon) < -180:
                            windLon += 360
                        allLonWs.append(windLon)
                        
    allIncDiffs = np.array(allIncDiffs)
    allLonDiffs = np.array(allLonDiffs)    
    
    # Split by model comparisons
    thetaW  = []
    thetaUW = []
    thetaA  = []
    thetaUA = []
    phiW    = []
    phiUW   = []
    phiA    = []
    phiUA   = []
    for i in range(int(data[-1,0])):
        j = i + 1
        idxW = np.where((data[:,0].astype(int) == j) & (data[:,3] == 'Wind'))[0]
        idxA = np.where((data[:,0].astype(int) == j) & (data[:,3] == 'ACE'))[0]
        if len(np.unique(data[idxW,4].astype(float))) > 1:
            # Inclination (the easy one)
            myIncs = data[idxW,4].astype(float)
            thetaW.append(np.mean(myIncs))
            thetaUW.append(np.std(myIncs))
            
            myLons = data[idxW,5].astype(float)
            if (np.max(myLons)-np.min(myLons)) > 180.:
                myLons[np.where(myLons > 180)] -= 360.
            phiW.append(np.mean(myLons))
            phiUW.append(np.std(myLons))
            
        if len(np.unique(data[idxA,4].astype(float))) > 1:
            # Inclination (the easy one)
            myIncs = data[idxA,4].astype(float)
            thetaA.append(np.mean(myIncs))
            thetaUA.append(np.std(myIncs))
            
            myLons = data[idxA,5].astype(float)
            if (np.max(myLons)-np.min(myLons)) > 180.:
                myLons[np.where(myLons > 180)] -= 360.
            phiA.append(np.mean(myLons))
            phiUA.append(np.std(myLons))
            
    thetaW  = np.array(thetaW)
    thetaUW = np.array(thetaUW)
    thetaA  = np.array(thetaA)
    thetaUA = np.array(thetaUA)
    phiW    = np.array(phiW)
    phiUW   = np.array(phiUW)
    phiA    = np.array(phiA)
    phiUA   = np.array(phiUA)
    
    # Just LLAMA
    data = np.genfromtxt('FRfits.dat', dtype=str)
    satTags = ['ACE','Wind']
    
    thetas = [[[] for i in range(4)],[[] for i in range(4)]]
    phis = [[[] for i in range(4)], [[] for i in range(4)]]
    
    dTllama = []
    dPllama = []
    for k in range(2):
        satTag = satTags[k]
        for i in range(int(data[-1,0])):
            j = i+1
            idx = np.where(data[:,0].astype(int) == j)
            subIdx = np.where((data[:,0].astype(int) == j) & (data[:,3]==satTag) )[0]
            names = np.unique(data[subIdx,2])
            isLlama = ['LLAMA' in name for name in names]
            llamaNames = names[isLlama]
            # Find the multi-Llama cases
            if (len(llamaNames) > 1) & ('LLAMAICEA' in llamaNames):
                thisCase = data[subIdx,:]
                # LA
                thisIdx = np.where((thisCase[:,2] == llamaNames[0]) & (thisCase[:,3]==satTag))[0]
                thetas[k][0].append(thisCase[thisIdx[0],4])
                phis[k][0].append(thisCase[thisIdx[0],5])
                # LB
                if 'LLAMAICEB' in llamaNames:
                    thisIdx = np.where((thisCase[:,2] == 'LLAMAICEB') & (thisCase[:,3]==satTag))[0]
                    thetas[k][1].append(thisCase[thisIdx[0],4])
                    phis[k][1].append(thisCase[thisIdx[0],5])
                else:
                    thetas[k][1].append(None)
                    phis[k][1].append(None)
                # LC
                if 'LLAMAICEC' in llamaNames:
                    thisIdx = np.where((thisCase[:,2] == 'LLAMAICEC') & (thisCase[:,3]==satTag))[0]
                    thetas[k][2].append(thisCase[thisIdx[0],4])
                    phis[k][2].append(thisCase[thisIdx[0],5])
                else:
                    thetas[k][2].append(None)
                    phis[k][2].append(None)
            
                # LD
                if 'LLAMAICED' in llamaNames:
                    thisIdx = np.where((thisCase[:,2] == 'LLAMAICED') & (thisCase[:,3]==satTag))[0]
                    thetas[k][3].append(thisCase[thisIdx[0],4])
                    phis[k][3].append(thisCase[thisIdx[0],5])
                else:
                    thetas[k][3].append(None)
                    phis[k][3].append(None)
    
    
        for i in range(4):
            thetas[k][i] = np.array(thetas[k][i]).astype(float)
            phis[k][i] = np.array(phis[k][i]).astype(float)
    
        diffsTheta = [[] for i in range(3)]
        diffsPhi = [[] for i in range(3)]
        # Clean up angles        
        for i in range(len(phis[k][0])):
            for jj in range(3):
                j = jj+1
                if phis[k][j][i] - phis[k][0][i]  > 180:
                    phis[k][j][i] -=360
                elif phis[k][j][i] - phis[k][0][i]  < -180:
                    phis[k][j][i] += 360
                #if j == 3:
                #    print (phis[j][i], phis[0][i])
            
                if np.isfinite(thetas[k][j][i]):
                    dTheta = np.abs(thetas[k][j][i] - thetas[k][0][i])
                    if dTheta > 90:
                        dTheta = 180 - dTheta
                    diffsTheta[jj].append(dTheta)
            
                if np.isfinite(phis[k][j][i]):
                    dPhi = np.abs(phis[k][j][i] - phis[k][0][i])
                    if dPhi > 180:
                        dPhi = 360 - dPhi
                    diffsPhi[jj].append(dPhi)
                
                if phis[k][j][i] < -80:
                    phis[k][j][i] += 360
                if phis[k][j][i] > 400:
                    phis[k][j][i] -=360
        dTllama.append(diffsTheta)
        dPllama.append(diffsPhi)
    
    
    
    fig, ax = plt.subplots(2,2, figsize=(7,7))
    ax = [ax[0,0], ax[0,1], ax[1,0], ax[1,1]]
    
    bmaxes = [90, 180]
    bins1 = np.linspace(0,bmaxes[0], 18)
    bins2 = np.linspace(0,bmaxes[1], 18)
    
    myParams = [[allIncDiffs, thetaUW, dTllama[1]], [allLonDiffs, phiUW, dPllama[1]], [allIncDiffs, thetaUA, dTllama[0]], [allLonDiffs, phiUA, dPllama[0]]]
    myBins = [bins1, bins2, bins1, bins2]
    labels = ['ACE/Wind', 'All Bounds', 'LLAMAICEB', 'LLAMAICEC', 'LLAMAICED']
    xlabs = ['$\\sigma_{\\theta,W}$', '$\\sigma_{\\phi,W}$', '$\\sigma_{\\theta,A}$', '$\\sigma_{\\phi,A}$']
    
    for j in range(4):
        hist, edges = np.histogram(myParams[j][0]/np.sqrt(2), bins=myBins[j], density=True )
        cents = 0.5*(edges[1:] + edges[:-1])
        myMax = np.max(hist)
        if j == 0:
            ax[j].plot(cents, hist/myMax, 'k', marker='o', label=labels[0])
        else:
            ax[j].plot(cents, hist/myMax, 'k', marker='o')

        hist, edges = np.histogram(myParams[j][1], bins=myBins[j], density=True )
        cents = 0.5*(edges[1:] + edges[:-1])
        if j == 0:
            ax[j].plot(cents, hist/myMax, 'gray', marker='o', label = labels[1] )
        else:
            ax[j].plot(cents, hist/myMax, 'gray', marker='o')
        Lnames = ['LLAMAICEB', 'LLAMAICEC', 'LLAMAICED']
        for i in range(3):
            hist, edges = np.histogram(myParams[j][2][i]/np.sqrt(2), bins=myBins[j], density=True )
            cents = 0.5*(edges[1:] + edges[:-1])
            if j == 0:
                ax[j].plot(cents, hist/myMax, catCols[Lnames[i]], marker='o', label = labels[2+i])
            else:
                ax[j].plot(cents, hist/myMax, catCols[Lnames[i]], marker='o')
        print (j, np.median(myParams[j][0]/np.sqrt(2)), np.median(myParams[j][1]))
        ax[j].set_xlim([0,myBins[j][-1]])
        ax[j].set_ylim([0,1])
        ax[j].set_ylabel('ScPDF')
        ax[j].set_xlabel(xlabs[j])
    #labels = [h.get_label() for h in handles]
    fig.legend( loc='upper center', ncol=3, fontsize=11)
    
    fig.subplots_adjust(hspace=0.3, wspace=0.3,  top=0.9, right=0.97, left=0.15, bottom=0.1)  
    plt.savefig('paperFigs/fitSPDFs.png')
    
#timelinePlot()
compAW()
#ratioPlot()
#catHistos()
#paramHistos()
#uncHistos()
#fitStats()
#compFitAW()
#compLLAMA()
#compRecons()
#compLLAMArecons()
#fitPDFs()