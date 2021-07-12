# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 18:42:19 2020

@author: ryanh
"""
# data managment and math functions
import pandas as pd
import numpy as np


import neuroseries as nts

# plotting
from matplotlib import pyplot as plt

# scipy
import scipy.io
import scipy.signal
from scipy import stats
from scipy.signal import hilbert,find_peaks
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
from ripple_detection import Karlsson_ripple_detector, filter_ripple_band

from ripple_detection.core import gaussian_smooth, get_envelope


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
    peak_time = []

    for ripple in ripple_times.itertuples():
        idx = np.logical_and(ts >= ripple.start_time, ts <= ripple.end_time)
        
        smooth_envelope = gaussian_smooth(get_envelope(filtered_lfps[idx,:]),0.004,fs)
        peaks = np.max(smooth_envelope,axis = 0)
        peak_idx = np.argmax(peaks)
        peak_time.append(ts[idx][np.argmax(smooth_envelope,axis=0)[peak_idx]])

        peak_amplitude.append(peaks[peak_idx])
        channel.append(peak_idx)
        
    ripple_times['peak_time'] = peak_time
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
        peakIx = scipy.signal.find_peaks(x = -rip, distance = peak_dist//(1/fs), threshold=0.0)
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

def emg_filter(session,ripple_times,shank,emg_thres=0.85):
    parts = session.split('\\')
    f = h5py.File(os.path.join(parts[0],parts[1],parts[2]) + '/EMG_from_LFP/' +
                  session.split('\\')[-1].split('.mat')[0] + '_emg.mat','r')
    emg = f['data'][0]
    emg_ts = f['timestamps'][0]
    max_emg=[]
    for ripple in ripple_times.itertuples():
        idx = np.logical_and(emg_ts >= ripple.start_time,
                             emg_ts <= ripple.end_time)
        if np.sum(idx) > 0:
            max_emg.append(np.max(emg[idx]))
        else:
            max_emg.append(1)
    
    ripple_times['max_emg'] = max_emg 
    
    if len(shank) > 8:
        ripple_times[np.array(max_emg) < emg_thres]
        
    return ripple_times

def make_Epochs(start, end):
    #Function to make an nts.IntervalSet dataframe with starting and ending epochs
    #Firstly, check whether both the lists are of same size or not
    if not (len(start) == len(end)):
        print("Start and End array lists are not of same dimension. Epochs IntervalSet can't be developed.")
        sys.exit()
    else:
        nts_array = []
        for i in range(len(start)):
            nts_array.append(nts.IntervalSet(start[i], end[i]))
        print(nts_array)
        return nts_array
    
def writeNeuroscopeEvents(path, ep, name):
    f = open(path, 'w')
    for i in range(len(ep)):
        f.writelines(str(ep.as_units('ms').iloc[i]['start']) + " "+name+" start "+ str(1)+"\n")
        #f.writelines(str(ep.as_units('ms').iloc[i]['peak']) + " "+name+" start "+ str(1)+"\n")
        f.writelines(str(ep.as_units('ms').iloc[i]['end']) + " "+name+" end "+ str(1)+"\n")
    f.close()
    return

def save_ripples(ripple_times,path):
    rpt_ep = nts.IntervalSet(np.array(ripple_times.start_time),
                             np.array(ripple_times.end_time),time_units = 's')

    writeNeuroscopeEvents(path + "\Swr_Ripple.evt.rip", rpt_ep, "SWR Ripple event")
     
def clipped(x, axis=1):
       x_diff = np.diff(x,axis=1)
       return np.sum(x_diff==0,axis=1) / x_diff.shape[1]
    
def clip_filter(ripple_times,ripple_maps,clip_thres=0.05):
    
    ripple_times['clipped'] = clipped(ripple_maps['ripple_map'])
    idx = ripple_times.clipped < clip_thres
    
    for key in ripple_maps.keys():
        ripple_maps[key] = ripple_maps[key][idx]
        
    ripple_times = ripple_times[idx]
    
    ripple_times= ripple_times.reset_index()
    ripple_times['ripple_number'] = np.arange(0,len(ripple_times),1)
    
    return ripple_times,ripple_maps
    
def filter_high_amp(ripple_times,ripple_maps,amp_thres=25):
    
    idx = ripple_times.peak_amplitude < amp_thres
    
    for key in ripple_maps.keys():
        ripple_maps[key] = ripple_maps[key][idx]
    
    ripple_times = ripple_times[idx]
    
    ripple_times= ripple_times.reset_index()
    ripple_times['ripple_number'] = np.arange(0,len(ripple_times),1)
    
    ripple_times = ripple_times.drop(columns=['index'])
    
    return ripple_times,ripple_maps

def filter_single_peaks(ripple_times,ripple_maps,peak_thres=0.30):
    peaks = []
    for x in ripple_maps['ripple_map']:
        # region around peak
        x = x[(len(x)//2 - 20) : (len(x)//2 + 20)]
        # center
        x = x - np.mean(x)
        # flip to greater mag 
        if np.abs(np.min(x)) > np.abs(np.max(x)):
            x = -x
        
        peak, _ = find_peaks(x,height=np.max(x)*peak_thres)
        peaks.append(len(peak))
    idx = np.array(peaks) > 1  
    
    for key in ripple_maps.keys():
        ripple_maps[key] = ripple_maps[key][idx]
    
    ripple_times = ripple_times[idx]
    
    ripple_times= ripple_times.reset_index()
    ripple_times['ripple_number'] = np.arange(0,len(ripple_times),1)
    
    ripple_times = ripple_times.drop(columns=['index'])
    return ripple_times,ripple_maps 
    

def get_good_channels(shank):
    #extract values from dictionary
    an_array = np.array(list(shank.values()),dtype=object)
    
    #loop through array to pull out individual channel        
    good_ch = []
    for i in range(len(an_array)):
        for x in range(len(an_array[i])):
            good_ch.append(an_array[i][x])
        
    return good_ch
    
def run_all(session):
    
    # get data session path from mat file
    path = get_session_path(session)
    
    # load position data from .mat file
    df = load_position(session)
    
    # load xml which has channel & fs info
    channels,fs,shank = loadXML(path)
    
    # get good channels
    good_ch = get_good_channels(shank)
    
    # load .lfp
    # lfp, ts = load_lfp(glob.glob(path +'\*.lfp')[0],channels,fs)
    lfp,ts = loadLFP(glob.glob(path +'\*.lfp')[0], n_channels=channels,
                     channel=good_ch, frequency=fs,
                     precision='int16')
    
    # interp speed of the animal
    speed = np.interp(ts,df.ts,df.speed)
    speed[np.isnan(speed)] = 0
    
    # detect ripples
    print('detecting ripples')
    ripple_times = Karlsson_ripple_detector(ts, lfp, speed, fs)
    
    # restrict ripples to < 200 ms
    ripple_times['ripple_duration'] = ripple_times.end_time - ripple_times.start_time
    ripple_times = ripple_times[ripple_times.ripple_duration < 0.200]
        
    # check against emg (< 0.85)
    ripple_times = emg_filter(session,ripple_times,shank)
        
    # get filtered signal
    print('filtering signal')
    LFPs = lfp
    filtered_lfps = np.stack([filter_ripple_band(lfp, fs) for lfp in LFPs.T])
    filtered_lfps = filtered_lfps.T
    
    # add ripple channel and peak amp
    print('getting ripple channel')
    ripple_times = get_ripple_channel(ripple_times,
                                      stats.zscore(filtered_lfps,axis=0),
                                      ts,fs)
    
    # get instant phase, amp, and freq
    print('get instant phase, amp, and freq')
    phase,amp,freq = get_phase_amp_freq(filtered_lfps,fs)
    
    # get ripple_map
    print('getting ripple maps')    
    ripple_maps = get_ripple_maps(ripple_times,ts,lfp,filtered_lfps,phase,amp,freq,fs)
    
    # get ripple frequency
    print('getting ripple frequency')    
    ripple_times['peak_freq'] = [map[len(map)//2] for map in ripple_maps['freq_map']]
    
    # filter out cliped signal
    ripple_times,ripple_maps = clip_filter(ripple_times,ripple_maps)
    
    # filter out very high amplitude ripples
    ripple_times,ripple_maps = filter_high_amp(ripple_times,ripple_maps)

    # find ripples with a single large jump
    ripple_times,ripple_maps = filter_single_peaks(ripple_times,ripple_maps)

    # save ripples for neuroscope inspection
    save_ripples(ripple_times,path)
    
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

# find HPC sessions
df_sessions = pd.read_csv('D:/ryanh/github\harvey_et_al_2020/Rdata_pae_track_cylinder_all_cells.csv')
sessions = pd.unique(df_sessions.session)
sessions = data_path+sessions

parallel = 0
# sessions.reverse()

if parallel==1:
    num_cores = multiprocessing.cpu_count()                             
    processed_list = Parallel(n_jobs=num_cores)(delayed(main_loop)(session,data_path,save_path) for session in sessions)
else:
    for session in sessions:
        sys.stdout.write('\rcurrent cell: %s' %(session))
        sys.stdout.flush()
        print(session)
        main_loop(session,data_path,save_path)
    

