import numpy as np
import numpy.matlib as npmatlib
import matplotlib.pyplot as plt
from scipy.fft import fftshift
from scipy.fft import fft
from scipy.fft import fftfreq
from scipy import signal
import math
import subprocess
from shutil import which
from scipy.signal import ShortTimeFFT
from scipy.signal.windows import hann
from scipy.signal.windows import boxcar
import os
from shutil import rmtree




def get_channel_psd(frequency_MHz = 150, duration=2, gain=21, plotLogic=False):
    complexData = get_airspy_raw_data(frequency_MHz, duration, gain)
    Fs = 3750;
    decimateStages = [8,5,5,2,2]
    dt = 1/Fs; 
    tWind = 0.02;
    nWind = math.floor(tWind * Fs);
    nOvlp = math.floor(tWind * Fs * 1/2);

    for stage in decimateStages:
        complexData = signal.decimate(complexData, stage)

    f, Pxx = signal.welch(complexData, Fs, window='hann', nperseg=nWind, noverlap=nOvlp,scaling='density',return_onesided=False,detrend=False)

    if plotLogic:
        figPSD, axPSD = plt.subplots(nrows = 1, ncols = 1, layout='constrained')
        nSamps = np.size(complexData)
        window = hann(nWind);
        scale  = 1/Fs * 1/np.sum(window**2)

        STF    = ShortTimeFFT(window,  hop = nWind-nOvlp, fs=Fs, fft_mode='twosided' )
        PxxMat = scale * np.abs(fftshift(STF.stft(complexData), axes=0))**2
        tSpec  = STF.t(nSamps)
        fSpec  = fftshift(STF.f)
        axPSD.pcolormesh(tSpec, fSpec, 10*np.log10(PxxMat), shading = 'gouraud' )

        figPSD2, axPSD2 = plt.subplots(nrows = 1, ncols = 1, layout='constrained')
        axPSD2.plot(fSpec, 10*np.log10(np.mean(PxxMat,axis=1)))
        axPSD2.plot(fftshift(f), fftshift(10*np.log10(Pxx))) 
        plt.show()
    
    return f, Pxx


def make_airspy_test_file(frequency_MHz=150, duration=2, gain=21, folderName=''):
    #Leave foldername empty if you want it in the current directory
    Fs = 3000000
    samps = math.ceil(Fs * duration)
    airspyrxInstallLoc = which('airspy_rx')
    #fileName = str(frequency_MHz)+ '_MHz_' + str(gain) + '_gain.bin'
    fileName =  f"{frequency_MHz:.6f}"+ '_MHz_' + str(gain) + '_gain.bin' #format string to always have 6 decimal places
    # folderName = 'RawData'
    fullName = os.path.join(folderName, fileName)
    airspyCommand = airspyrxInstallLoc + ' -f ' + str(frequency_MHz) + ' -a 3000000 -t 0 -h ' + str(gain) + ' -n ' + str(samps) + ' -r ' + fullName
    p1 = subprocess.run(airspyCommand,shell=True,capture_output=True,timeout=duration+5)
    return p1.returncode, fileName


def get_airspy_raw_data(frequency_MHz = 150, duration=2, gain=21, dataFolderName=''):
    leadTime = 0.1
    Fs       = 3000000
    status, fileName = make_airspy_test_file(frequency_MHz, duration+leadTime, gain, dataFolderName)
    fullName = os.path.join(dataFolderName, fileName)
    if status != 0 :
        raise ChildProcessError('Airspy data collection failed. ')
    #fileName    = "test.bin"
    dataRaw     = np.fromfile(fullName, dtype=np.float32 )
    complexData = np.array(dataRaw[:-1:2] + 1j*dataRaw[1::2])
    nSamps2Cut  = int(np.ceil( leadTime * Fs ) )
    complexData = complexData[nSamps2Cut::]
    #nSampsRaw   = np.size(complexData)
    return complexData
    

def get_new_airspy_psd(frequency_MHz=150, duration=2, gain=21, dataFolderName=''):
    complexData = get_airspy_raw_data(frequency_MHz, duration, gain, dataFolderName)
    FsRaw = 3000000;
    tWind = 0.02;
    nWind = math.floor(tWind * FsRaw);
    nOvlp = math.floor(tWind * FsRaw * 1/2);
    f, Pxx = signal.welch(complexData, FsRaw, window='hann', nperseg=nWind, noverlap=nOvlp,scaling='density',return_onesided=False,detrend=False)
    return f, Pxx

def get_existing_airspy_psd(fileName, Fs):
    dataRaw     = np.fromfile(fileName, dtype=np.float32 )
    complexData = np.array(dataRaw[:-1:2] + 1j*dataRaw[1::2])
    tWind = 0.02;
    nWind = math.floor(tWind * Fs);
    nOvlp = math.floor(tWind * Fs * 1/2);
    f, Pxx = signal.welch(complexData, Fs, window='hann', nperseg=nWind, noverlap=nOvlp,scaling='density',return_onesided=False)
    # if plotLogic:
    #     figPSD, axPSD = plt.subplots(nrows = 1, ncols = 1, layout='constrained')
    #     nSamps = np.size(complexData)
    #     window = hann(nWind);
    #     STF    = ShortTimeFFT(window,  hop = nWind-nOvlp, fs=Fs, fft_mode='twosided' )
    #     Pxx    = (1/Fs)**2/tWind * np.abs(fftshift(STF.stft(complexData), axes=0))**2
    #     t      = STF.t(nSamps)
    #     f      = fftshift(STF.f)
    #     axPSD.pcolormesh(t, f, 10*np.log10(Pxx), shading = 'gouraud' )
    return f, Pxx


def get_multiple_airspy_psd(centerFreqs_MHz = 150, duration=2, gain=21, plotLogic=False, dataFolderName=''):
    firstRun = True
    if dataFolderName:#check to see if the folder name was supplied. if not, we don't need a directory and the files will be written to the current working directory
        try:
            os.mkdir(dataFolderName)
        except:
            #Directory already exists - so clear it out
            print('Data folder already exists. Clearing contents. ')
            clear_directory(dataFolderName)
    

    if plotLogic:
        figFFT, axFFT = plt.subplots(nrows = 1, ncols = 1, layout='constrained')
    for indx, freq in enumerate(centerFreqs_MHz): #https://stackoverflow.com/questions/522563/how-to-access-the-index-value-in-a-for-loop
        print("Collecting frequency: " + str(freq) + " MHz")
        f, Pxx = get_new_airspy_psd(freq, duration, gain, dataFolderName)
        if firstRun:
            PxxArray = np.zeros((np.size(Pxx), np.size(centerFreqs_MHz)))
            firstRun = False
        PxxArray[:, indx] = Pxx
        if plotLogic:
            axFFT.plot(fftshift(f + freq*10**6), fftshift(10*np.log10(Pxx)))
    return f, PxxArray


#Written by chatGPT with query:
#write a python function that takes a directory as an input string and then deletes the contents of that directory
def clear_directory(directory_path):
    """
    Deletes all the contents of the given directory.

    Args:
    directory_path (str): Path to the directory to clear.
    """
    if not os.path.exists(directory_path):
        print(f"Directory '{directory_path}' does not exist.")
        return

    if not os.path.isdir(directory_path):
        print(f"'{directory_path}' is not a directory.")
        return

    # Iterate through all the files and folders in the directory
    for item_name in os.listdir(directory_path):
        item_path = os.path.join(directory_path, item_name)

        try:
            # If the item is a file, remove it
            if os.path.isfile(item_path) or os.path.islink(item_path):
                os.unlink(item_path)
            # If the item is a directory, remove it and all its contents
            elif os.path.isdir(item_path):
                rmtree(item_path)
        except Exception as e:
            print(f"Failed to delete {item_path}. Reason: {e}")

# Example usage
# clear_directory("/path/to/your/directory")


################################################
# OLD DEVELOPMENT CODE
################################################

# t = time.time()
# #plt.ioff()  # Use non-interactive mode.
# testDescription = "In_Office_Testing"
# fRange_MHz      = [146, 152]##146.539#[146, 148]
# duration        = 1
# freqStep_MHz    = 0.375

# # figFFT, axFFT = plt.subplots(nrows = 1, ncols = 1, layout='constrained')
# # fChan, PxxChan = get_channel_psd(148, duration, 21, False)
# # axFFT.plot(fftshift(fChan + 148.00*10**6), fftshift(10*np.log10(PxxChan)))

# # complexData = get_airspy_raw_data(148, 2, 21)
# # figComplex, axComplex = plt.subplots(nrows = 1, ncols = 1, layout='constrained')
# # axComplex.plot(1/3000000*np.arange(np.size(complexData)), np.real(complexData))

# # for gain in gainSettings[::4]:
# #     fChan, PxxChan = get_channel_psd(148.00, duration, gain, False)
# #     axFFT.plot(fftshift(fChan + 148.00*10**6), fftshift(10*np.log10(PxxChan)))


# make_airspy_test_file(150, 30, 21)#Run radio for 30s to warm up 
# nFreqSteps   = 1 + (np.max(fRange_MHz) - np.min(fRange_MHz)) / (freqStep_MHz)
# centerFreqs = np.min(fRange_MHz) + freqStep_MHz * np.arange(nFreqSteps)
# # gainSettings = np.arange(22)
# gainSettings = [10, 12, 14, 16, 18, 20, 21]

# figFFTMulti, axFFTMulti = plt.subplots(nrows = 1, ncols = 1, layout='constrained')

# for gain in gainSettings:
#     f, PxxArray = get_multiple_airspy_psd(centerFreqs, duration, gain)
#     fArrayCenters, fArrayDiffs = np.meshgrid(centerFreqs, f) 
#     fArray   = np.round((fArrayCenters * 10**6 + fArrayDiffs) / 10**6, 6)
#     fShift   = fftshift(f)
#     fSelectMask = np.logical_and(fShift>=-freqStep_MHz/2*10**6, fShift<freqStep_MHz/2*10**6)
#     PxxArray = fftshift(PxxArray,axes=0)
    
#     np.savetxt("Pxx_W_Hz_"+ testDescription + '_Gain=' + str(gain) + ".csv",PxxArray, delimiter=",")
#     for indx in np.arange(np.size(centerFreqs)):
#         #axFFT.plot(fftshift(f + centerFreqs[indx]*10**6), fftshift(10*np.log10(PxxArray[:,indx])))
#         axFFTMulti.plot(fShift[fSelectMask] + centerFreqs[indx]*10**6, (10*np.log10(PxxArray[fSelectMask,indx])))

#     plt.show()    


# np.savetxt("F_MHz_"+ testDescription +".csv",fArray, delimiter=",",fmt='%3.6f')

# elapsed = time.time() - t

# print('Execution Time: ' + str(elapsed))


    


# fDogTrack, PxxDogTrack = get_existing_airspy_psd('/Users/mshafer/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/FLIGHT_TESTING_DATA/2023-12-14-Zimbabwe Thursday/09-Airspy capture data f146539000-Peters antenna/2023-12-14-Dog_tracking_at_16_16-f146539000.bin', 3000000)
# axFFT.plot(fftshift(fDogTrack + 146.539*10**6), fftshift(10*np.log10(PxxDogTrack)))
# axFFT.plot(fftshift(fDogTrack + 146.539*10**6), fftshift(10*np.log10(PxxDogTrack)))
# fDogTrackFlight, PxxDogTrackFlight = get_existing_airspy_psd('/Users/mshafer/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/FLIGHT_TESTING_DATA/2023-12-14-Zimbabwe Thursday/08-Flight 4 Dog Tracking/RPI_LOGS_23-12-14_22_26_24/data_record.2.1.bin', 3750)
# axFFT.plot(fftshift(fDogTrackFlight + 146.539*10**6), fftshift(10*np.log10(PxxDogTrackFlight)))


# plt.show
# 
# fftFreqs, psd = get_airspy_noise(146.10, True, False)
# axFFT[1].plot(fftFreqs, 10*np.log10(psd))



#plt.show
   


# from scipy.signal import ShortTimeFFT
# from scipy.signal.windows import hann

# decimateStages = [8,5,5,2,2]
# for stage in decimateStages:
#     complexData = signal.decimate(complexData, stage)

# nSamps = np.size(complexData)

# Fs = 3750;
# dt = 1/Fs; 
# tWind = 0.02;
# nWind = math.floor(tWind * Fs);
# nOvlp = math.floor(tWind * Fs * 1/2);
# time = np.arange(nSamps) * dt

# fig, ax = plt.subplots(nrows = 2, ncols = 1, layout='constrained')

# ax[0].plot(time, complexData.real)

# window = hann(nWind);

# STF = ShortTimeFFT(window,  hop = nWind-nOvlp, fs=Fs, fft_mode='twosided' )
# Sx = fftshift(STF.stft(complexData), axes=0)
# t = STF.t(nSamps)
# f = fftshift(STF.f)

# ax[1].pcolormesh(t, f, 10*np.log10(np.abs(Sx)**2), shading = 'gouraud' )

# meanSpectralDensity = 1/tWind * np.mean(np.abs(Sx)**2,axis=1)

# fig2, ax2 = plt.subplots(nrows = 1, ncols = 1, layout='constrained')
# ax2.plot(f, 10*np.log10(meanSpectralDensity))

# coreThreshold = 0.8
# fMin = np.min(f)
# fMax = np.max(f)
# coreMask = (f>=coreThreshold*fMin) & (f<=coreThreshold * fMax)
# corePSD_dB = 10*np.log10(np.mean(meanSpectralDensity[coreMask]))
# print('Core mean noise level is: ' + str(corePSD_dB) + ' dB')

# plt.show()




# plt.show()

# f, t, Zxx = signal.stft(complexData, Fs, window = 'hann', nperseg=nWind, noverlap=nOvlp , return_onesided = False)
# f = fftshift(f);
# Zxx = fftshift(Zxx, axes=0)
# np.savetxt("Zxx.csv", Zxx, delimiter = ",")
#ax[1].pcolormesh(t, f, 10*np.log10(np.abs(Zxx)**2), shading = 'gouraud' )



# def moving_average(a, n=3):
#     #Found from https://stackoverflow.com/questions/14313510/how-to-calculate-rolling-moving-average-using-python-numpy-scipy
#     #Modified to provide averages for the first n-1 samples as well with smaller window
#     n = int(n)
#     ret = np.cumsum(a, dtype=float)
#     ret[n:] = ret[n:] - ret[:-n]
#     #return ret[n - 1:]/n
#     ret[n - 1:] = ret[n - 1:]/n
#     ret[:n - 1] = ret[:n - 1]/(1 + np.arange(n-1))
#     return ret



    # fftData     = fftshift(fft(complexData))
    # fftFreqs    = fftshift(fftfreq(nSampsRaw, 1/FsRaw))
    # tAllData    = 1/FsRaw * nSampsRaw 
    # dt          = 1/FsRaw; 
    # psdAll      = dt**2/tAllData * np.abs(fftData)**2
    # freqWind    = 3750
    # nMovMean    = int(FsRaw/freqWind)

    # # figFFT, axFFT = plt.subplots(nrows = 2, ncols = 1, layout='constrained')
    # # axFFT[0].plot(fftFreqs, 10*np.log10(psdAll))
    # # axFFT[1].plot(fftFreqs, 10*np.log10(psdMovMean))
    # # plt.show
    # if movMean:
    #     psdAll  = moving_average(psdAll, nMovMean)

    # if resample:
    #     psdAll, fftFreqs   = signal.resample(psdAll,  nMovMean, fftFreqs)

    # return fftFreqs, psdAll, f, Pxx
