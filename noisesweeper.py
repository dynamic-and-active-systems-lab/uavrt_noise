import noisereader as nr
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.fft import fftshift
from os import path
import sys


def noisesweeper(testDescription):
    t = time.time()

    # print('Running radio for 10s to warm up...')
    # nr.make_airspy_test_file(150, 10, 21)#Run radio to warm up 
    # print('complete.')

    fRange_MHz      = [146, 154]#[146, 152]##146.539#[146, 148]
    duration        = 2
    freqStep_MHz    = 0.375
    gainSettings    = [21]


    nFreqSteps   = 1 + (np.max(fRange_MHz) - np.min(fRange_MHz)) / (freqStep_MHz)
    centerFreqs = np.min(fRange_MHz) + freqStep_MHz * np.arange(nFreqSteps)
    

    figFFTMulti, axFFTMulti = plt.subplots(nrows = 1, ncols = 1, layout='constrained')

    for gain in gainSettings:
        print('Recording gain setting: ' + str(gain))
        folderName = 'Data-' + testDescription
        print('The folder name is: '+ folderName)
        f, PxxArray = nr.get_multiple_airspy_psd(centerFreqs, duration, gain, False, folderName)
        fArrayCenters, fArrayDiffs = np.meshgrid(centerFreqs, f) 
        fArray      = np.round((fArrayCenters * 10**6 + fArrayDiffs) / 10**6, 6)
        #fArrayShift = fftshift(fArray,axes=0)
        fShift      = fftshift(f)
        fSelectMask = np.logical_and(fShift>=-freqStep_MHz/2*10**6, fShift<freqStep_MHz/2*10**6)
        pxxFullName = path.join(folderName, "Pxx_W_Hz_"+ testDescription + '_Gain=' + str(gain) + ".csv")
        np.savetxt(pxxFullName, PxxArray, delimiter=",")

        PxxArray    = fftshift(PxxArray,axes=0)

        for indx in np.arange(np.size(centerFreqs)):
            #axFFT.plot(fftshift(f + centerFreqs[indx]*10**6), fftshift(10*np.log10(PxxArray[:,indx])))
            axFFTMulti.plot(fShift[fSelectMask] + centerFreqs[indx]*10**6, (10*np.log10(PxxArray[fSelectMask,indx])))
    
    fFullName = path.join(folderName, "F_MHz_"+ testDescription +".csv")
    np.savetxt(fFullName, fArray, delimiter=",",fmt='%3.6f')

    plt.show()    


    elapsed = time.time() - t

    print('Execution Time: ' + str(elapsed))

if __name__=="__main__":
    testDescription = str(sys.argv[1])
    noisesweeper(testDescription)