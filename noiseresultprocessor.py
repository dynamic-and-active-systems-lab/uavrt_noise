import numpy as np
from scipy.fft import fftshift
import matplotlib.pyplot as plt
import os.path


def moving_average(a, n=3):
    #Found from https://stackoverflow.com/questions/14313510/how-to-calculate-rolling-moving-average-using-python-numpy-scipy
    #Modified to provide averages for the first n-1 samples as well with smaller window
    n = int(n)
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    #return ret[n - 1:]/n
    ret[n - 1:] = ret[n - 1:]/n
    ret[:n - 1] = ret[:n - 1]/(1 + np.arange(n-1))
    return ret


######################################################
# MODIFY THE FOLLOWING LINES OF CODE FOR YOUR SPECIFIC TEST
gainLevels      = [21] #[10, 12, 14, 16, 18, 20, 21]
folderStr       = "" #sub directory where the data files live or the full path to the directory where they live
testFileStrList = ["Testing"]#["OutDoorMacRecAntOnLadder", "OutdoorTestRPIRecAntOnLadderWifiAway","OutdoorTestRPIRecAntOnLadderDroneOnPoweringRpi","OutdoorTestRPIRecAntOnLadderDroneWithAnt","OutdoorTestRPIRecDroneFlying"]
colorsForTest   = ["blue"] #['blue','green', 'black', '#00FF00','red']
descriptions    = ["A simple test"]#['Ant on Ladder, Mac On Ground','Ant on Ladder, RPI On Ground w/ Ext. Power ', 'Ant on Ladder, RPI On Ground w/ Drone Power ', 'Drone with Ant on Ladder', 'Drone flying']
######################################################

baseFileStr     = "Pxx_W_Hz_"
freqFileStr     = "F_MHz_"
gainFileStr     = "_Gain="
extFileStr      = ".csv"

nGainLevls = len(gainLevels)

fig, ax = plt.subplots(nrows = nGainLevls, ncols = 1, layout='constrained')
figMov, axMov = plt.subplots(nrows = nGainLevls, ncols = 1, layout='constrained')

testNum = 0
for testFileStr in testFileStrList:
    print("Processing " + testFileStr )
    
    gainNum = 0
    for gain in gainLevels:
        fileStr      = os.path.join(folderStr, baseFileStr + testFileStr + gainFileStr + str(gain) + extFileStr)
        PxxMat       = np.loadtxt(fileStr, delimiter=",")
        PxxMat       = fftshift(PxxMat,axes=0)
        fileStr      = os.path.join(folderStr, freqFileStr + testFileStr + extFileStr)
        freqMat      = np.loadtxt(fileStr, delimiter=",")
        centerFreqs  = freqMat[0,:]
        freqStep     = centerFreqs[1] - centerFreqs[0]
        freqMat      = fftshift(freqMat, axes=0)
        indxMask     = np.logical_and( freqMat[:,0] >= centerFreqs[0]-freqStep/2,  freqMat[:,0] < centerFreqs[0]+freqStep/2)
        freqMatCore  = freqMat[indxMask, :]
        PxxMatCore   = PxxMat[indxMask, :]
        nRows, nCols = np.shape(PxxMatCore)
        colNum = 0
        
        ax.plot(freqMatCore.flatten(order='F'), 10*np.log10(PxxMatCore.flatten(order='F')), color=colorsForTest[testNum], label=descriptions[testNum])
        axMov.plot(moving_average(freqMatCore.flatten(order='F'),75), 10*np.log10(moving_average(PxxMatCore.flatten(order='F'),75)), color=colorsForTest[testNum], label=descriptions[testNum])

        # for Pxxcolumn in PxxMatCore.T:
        #     print("Plotting column " + str(colNum))
        #     #ax[gainNum].plot(freqMatCore[:, colNum], 10*np.log10(Pxxcolumn))
        #     ax.plot(freqMatCore[:, colNum], 10*np.log10(Pxxcolumn), color=colorsForTest[testNum], label=descriptions[testNum])
        #     axMov.plot(moving_average(freqMatCore[:, colNum],75),10*np.log10(moving_average(Pxxcolumn,75)), color=colorsForTest[testNum], label=descriptions[testNum])
        #     colNum += 1
        gainNum += 1
    testNum += 1
   # plt.close(fig)

ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel('Power spectral density (W/Hz)')
ax.grid(True)
ax.legend()
axMov.set_xlabel('Frequency (MHz)')
axMov.set_ylabel('Power spectral density (W/Hz)')
axMov.grid(True)
axMov.legend()

fig.savefig(folderStr+".pdf", bbox_inches='tight')
figMov.savefig(folderStr + "_MovMean" + ".pdf", bbox_inches='tight')
plt.show()




        


