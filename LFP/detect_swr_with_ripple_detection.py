# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 18:42:19 2020

@author: ryanh
"""
# data managment and math functions
import pandas as pd
import numpy as np
import math

import time

# import neuroseries as nts

# plotting
from matplotlib import pyplot as plt

# scipy
import scipy.io
import scipy.signal
from scipy import stats
from scipy.signal import hilbert,detrend
from scipy.ndimage import gaussian_filter1d

# for loading files
import h5py
import sys,os
import glob
import pickle

# parallel processing
import multiprocessing
from joblib import Parallel, delayed

# ripple detector
from ripple_detection import Kay_ripple_detector, Karlsson_ripple_detector, filter_ripple_band

from ripple_detection.core import gaussian_smooth, get_envelope

# Continuous Wavelet Transform
import obspy
from obspy.imaging.cm import obspy_sequential
from obspy.signal.tf_misfit import cwt


def loadXML(path):
    """
    path should be the folder session containing the XML file
    Function returns :
        1. the number of channels
        2. the sampling frequency of the dat file or the eeg file depending of what is present in the folder
            eeg file first if both are present or both are absent
        3. the mappings shanks to channels as a dict
    Args:
        path : string
    Returns:
        int, int, dict
    """
    if not os.path.exists(path):
        print("The path "+path+" doesn't exist; Exiting ...")
        sys.exit()
    listdir = os.listdir(path)
    xmlfiles = [f for f in listdir if f.endswith('.xml')]
    if not len(xmlfiles):
        print("Folder contains no xml files; Exiting ...")
        sys.exit()
    new_path = os.path.join(path, xmlfiles[0])

    from xml.dom import minidom
    xmldoc = minidom.parse(new_path)
    nChannels = xmldoc.getElementsByTagName('acquisitionSystem')[0].getElementsByTagName('nChannels')[0].firstChild.data
    fs_dat = xmldoc.getElementsByTagName('acquisitionSystem')[0].getElementsByTagName('samplingRate')[0].firstChild.data
    fs = xmldoc.getElementsByTagName('fieldPotentials')[0].getElementsByTagName('lfpSamplingRate')[0].firstChild.data

    shank_to_channel = {}
    groups = xmldoc.getElementsByTagName('anatomicalDescription')[0].getElementsByTagName('channelGroups')[0].getElementsByTagName('group')
    for i in range(len(groups)):
        shank_to_channel[i] = np.sort([int(child.firstChild.data) for child in groups[i].getElementsByTagName('channel')])
    return int(nChannels), int(fs), shank_to_channel

def loadLFP(path, n_channels=90, channel=64, frequency=1250.0, precision='int16'):
    if type(channel) is not list:
        f = open(path, 'rb')
        startoffile = f.seek(0, 0)
        endoffile = f.seek(0, 2)
        bytes_size = 2
        n_samples = int((endoffile-startoffile)/n_channels/bytes_size)
        duration = n_samples/frequency
        interval = 1/frequency
        f.close()
        with open(path, 'rb') as f:
            data = np.fromfile(f, np.int16).reshape((n_samples, n_channels))[:,channel]
            timestep = np.arange(0, len(data))/frequency
            return data, timestep # nts.Tsd(timestep, data, time_units = 's')
        
    elif type(channel) is list:
        f = open(path, 'rb')
        startoffile = f.seek(0, 0)
        endoffile = f.seek(0, 2)
        bytes_size = 2

        n_samples = int((endoffile-startoffile)/n_channels/bytes_size)
        duration = n_samples/frequency
        f.close()
        with open(path, 'rb') as f:
            data = np.fromfile(f, np.int16).reshape((n_samples, n_channels))[:,channel]
            timestep = np.arange(0, len(data))/frequency
            return data,timestep # nts.TsdFrame(timestep, data, time_units = 's')
        
def load_lfp(file,channels,fs):
    LFP = []
    for i in range(channels):
        lfp = loadLFP(file, n_channels=channels, channel=i, frequency=fs)
        LFP.append(lfp)

    lfp = np.vstack(LFP)  
    ts = np.arange(0, lfp.shape[1])/fs
    return lfp.T, ts

def load_position(session):
    f = h5py.File(session,'r')
    # load frames [ts x y a s] 
    frames = np.transpose(np.array(f['frames']))
    return pd.DataFrame(frames,columns=['ts', 'x', 'y', 'hd', 'speed'])    

def get_ripple_channel(ripple_times,filtered_lfps,ts,fs):
    channel = []
    peak_amplitude = []

    for ripple in ripple_times.itertuples():
        idx = np.logical_and(ts >= ripple.start_time, ts <= ripple.end_time)
        
        smooth_envelope = gaussian_smooth(get_envelope(filtered_lfps[idx,:]),0.004,fs)
        peaks = np.max(smooth_envelope,axis = 0)
        peak_idx = np.argmax(peaks)

        peak_amplitude.append(peaks[peak_idx])
        channel.append(peak_idx)

    ripple_times['peak_channel'] = channel
    ripple_times['peak_amplitude'] = peak_amplitude
    return ripple_times

def get_phase_amp_freq(sig,fs):
    
    phas = []
    amp = []
    freq = []
    
    for signal in sig.T:
        analytic_signal = hilbert(signal)
        amplitude_envelope = np.abs(analytic_signal)
        phase = np.angle(analytic_signal)
        instantaneous_phase = np.unwrap(phase)
        instantaneous_frequency = gaussian_filter1d((np.diff(instantaneous_phase) / (2.0*np.pi) * fs),
                                                    0.004 * fs, truncate=8, axis=0,mode='constant')
        phas.append(phase)
        amp.append(amplitude_envelope)
        freq.append(instantaneous_frequency)
    
    phas = np.vstack(phas) 
    amp = np.vstack(amp)
    freq = np.vstack(freq) 

    return phas.T,amp.T,freq.T

def get_ripple_freq(ripple_times,freq,dt):
    peak_freq = []
    for ripple in ripple_times.itertuples():
        idx = np.logical_and(dt >= ripple.start_time, dt <= ripple.end_time)
        rip = freq[idx,ripple.peak_channel]
        peak_freq.append(rip[len(rip) // 2])
    ripple_times['peak_freq'] = peak_freq
    return ripple_times

def get_ripple_freq_peaks_method(ripple_times,filtered_lfps,ts,fs,peak_dist=0.0032):
    fqcy = np.zeros((len(ripple_times),1))
    i = 0
    for ripple in ripple_times.itertuples():
        idx = np.logical_and(ts >= ripple.start_time, ts <= ripple.end_time)
        rip = filtered_lfps[idx,ripple.peak_channel]
        # find peaks with a distance of 3.2 ms
        peakIx = scipy.signal.find_peaks(x = -rip, distance = peak_dist//(1/fs), threshold = 0.0)
        peakIx = peakIx[0]
        if (not (peakIx.size == 0)) and (peakIx.size != 1):
            fqcy[i] = fs/np.median(np.diff(peakIx))
        i += 1
    ripple_times['peak_freq'] = fqcy
    return ripple_times
    
def get_session_path(session):
    f = h5py.File(session,'r')
    return f['session_path'][()].tobytes()[::2].decode()

def get_ripple_maps(ripple_times,ts,lfp,filtered_lfps,phase,amp,freq,fs):
    
    # Initializing variables
    rip = np.zeros((len(ripple_times),151))
    rip_filt = np.zeros((len(ripple_times),151))
    rip_phase = np.zeros((len(ripple_times),151))
    rip_amp = np.zeros((len(ripple_times),151))
    rip_freq = np.zeros((len(ripple_times),151))

    # row index
    ind = np.arange(0,len(lfp),1) 

    i = 0
    for ripple in ripple_times.itertuples():
        # get ripple index
        idx = np.logical_and(ts >= ripple.start_time, ts <= ripple.end_time)
        # find peak of ripple using the smoothed filtered signal
        smooth_envelope = gaussian_smooth(get_envelope(filtered_lfps[idx,int(ripple.peak_channel)]),0.004,fs)
        rip_peak_idx = np.argmax(smooth_envelope)
        # find that peaks location in signal
        middle_idn = ind[idx][rip_peak_idx]
        # create expanded index
        idx = np.arange(middle_idn - 75,middle_idn + 76,1)
        
        # if ripple is the the very beginning or end of session
        if (middle_idn - 75 < 0) or (middle_idn + 76 > len(ind)):
            x = np.zeros(151)
            rip[i] = x
            rip_filt[i] = x
            rip_phase[i] = x
            rip_amp[i] = x
            rip_freq[i] = x
            print('ripple close to edge of session')
        else:
            # pull out expanded index
            rip[i] = lfp[idx,ripple.peak_channel]
            rip_filt[i] = filtered_lfps[idx,ripple.peak_channel]
            rip_phase[i] = phase[idx,ripple.peak_channel]
            rip_amp[i] = amp[idx,ripple.peak_channel]
            rip_freq[i] = freq[idx,ripple.peak_channel]

            i+=1
        
    ripple_maps = {"ripple_map": rip,
            "filtered_map":rip_filt,
            "phase_map":rip_phase,
           "amp_map":rip_amp,
           "freq_map":rip_freq}
    
    return ripple_maps

def normalize(list, range):
    l = np.array(list) 
    a = np.max(l)
    c = np.min(l)
    b = range[1]
    d = range[0]
    m = (b - d) / (a - c)
    pslope = (m * (l - c)) + d
    return pslope

def get_scalogram(sig,fs,padding=100,f_min=150,f_max=250,fig=1,ax=0):
    
    # sample difference
    dt = 1/fs
    # pad signal
    sig_padded = np.pad(sig, (padding, padding), 'linear_ramp')
    # get time stamps
    t = np.linspace(0, dt * len(sig), len(sig))
    # get scalogram
    scalogram = cwt(sig_padded, dt, 8, f_min, f_max)
    # delete padding
    scalogram = np.delete(scalogram, np.s_[1:padding+1], axis=1) 
    scalogram = np.delete(scalogram, np.s_[-(padding+1):-1], axis=1) 
    
    # plot figure
    if fig==1:
        if ax == 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        x, y = np.meshgrid(
            t,
            np.logspace(np.log10(f_min), np.log10(f_max), scalogram.shape[0]))

        ax.pcolormesh(x, y, np.abs(scalogram), cmap=obspy_sequential,shading='auto')
        ax.plot(t,normalize(sig,[f_min,f_max]),'k',linewidth=1)
        ax.set_ylabel("Frequency [Hz]")
        ax.set_ylim(f_min, f_max)
        # ax.set_yscale('log')
        if ax == 0:
            plt.show()
    
    return np.abs(scalogram)

def run_all(session):
    
    # get data session path from mat file
    path = get_session_path(session)
    
    # load position data from .mat file
    df = load_position(session)
    
    # load xml which has channel & fs info
    channels,fs,shank = loadXML(path)
    
    # load .lfp
    # lfp, ts = load_lfp(glob.glob(path +'\*.lfp')[0],channels,fs)
    lfp,ts = loadLFP(glob.glob(path +'\*.lfp')[0], n_channels=channels, channel=np.arange(0,channels,1), frequency=1250.0, precision='int16')
    
    # interp speed of the animal
    speed = np.interp(ts,df.ts,df.speed)
    speed[np.isnan(speed)] = 0
    
    # detect ripples
    print('detecting ripples')
    ripple_times = Karlsson_ripple_detector(ts, lfp, speed, fs,                       
                            speed_threshold=3.0, minimum_duration=0.02,
                            zscore_threshold=3.0, smoothing_sigma=0.004,
                            close_ripple_threshold=0.0)
    
    # restrict ripples to < 150 ms
    ripple_times['ripple_duration'] = ripple_times.end_time - ripple_times.start_time
    ripple_times = ripple_times[ripple_times.ripple_duration < 0.150]
    
    # get filtered signal
    print('filtering signal')
    LFPs = lfp
    filtered_lfps = np.stack([filter_ripple_band(lfp, fs) for lfp in LFPs.T])
    filtered_lfps = filtered_lfps.T
    
    # add ripple channel and peak amp
    print('getting ripple channel')
    ripple_times = get_ripple_channel(ripple_times,filtered_lfps,ts,fs)
    

    # get instant phase, amp, and freq
    print('get instant phase, amp, and freq')
    phase,amp,freq = get_phase_amp_freq(filtered_lfps,fs)
    dt = ts[:-1] + (1/fs)/2
    
    # get ripple frequency
    print('getting ripple frequency')
    ripple_times = get_ripple_freq(ripple_times,freq,dt)
    
    # get ripple_map
    print('getting ripple maps')    
    ripple_maps = get_ripple_maps(ripple_times,ts,lfp,filtered_lfps,phase,amp,freq,fs)

    return ripple_times,lfp,filtered_lfps,ts,ripple_maps

        
def main_loop(session,data_path,save_path):
    base = os.path.basename(session)
    os.path.splitext(base)
    save_file = save_path + os.path.splitext(base)[0] + '.pkl'
    
    # check if saved file exists
    if os.path.exists(save_file):
        return
        
    # detect ripples and calc some features
    ripple_times,lfp,filtered_lfps,ts,ripple_maps = run_all(session)   

    # save file
    with open(save_file, 'wb') as f:
        pickle.dump(ripple_times, f)
        pickle.dump(ripple_maps, f)


data_path = 'F:\\Projects\\PAE_PlaceCell\\ProcessedData\\'
save_path = "F:\\Projects\\PAE_PlaceCell\\swr_data\\"
sessions = glob.glob(data_path + '*.mat')
# sessions.reverse()

for session in sessions:
#     sys.stdout.write('\rcurrent cell: %s' %(session))
#     sys.stdout.flush()
    print(session)
    main_loop(session,data_path,save_path)
    
# num_cores = multiprocessing.cpu_count()                             
# processed_list = Parallel(n_jobs=num_cores)(delayed(main_loop)(session,data_path,save_path) for session in sessions)
