#%% This functions are custom made for preprocessing 
#   
from scipy.signal import medfilt,iirnotch,tf2zpk,zpk2sos,welch,windows
from custom_recording_extractors import load_notches_structs
import numpy as np
from spikeinterface.preprocessing.filter import  FilterRecording as filter
from spikeinterface import ChannelSliceRecording,ChannelsAggregationRecording 


#%%
def filter_freq_peaks(recording, dtype=np.float32, load_mat_notches=True, **kwargs):
    if load_mat_notches:
        # to load old mat file with notches
        notch_info = load_notches_structs(recording)
    else:
        # maybe somehting like ,**kwargs can be used to pass optional inputs about conpute notches(FC)
        notch_info = compute_notches_structs(recording, **kwargs)
    reclist = []

    for chIDs in  recording.get_channel_ids():
        subrec = ChannelSliceRecording(recording, [chIDs])
        if chIDs in notch_info:
            sos = notch_info[chIDs]['sos']
            sos[0,0:3] = sos[0,0:3] * notch_info[chIDs]['g']
            reclist.append(filter(subrec, coeff=sos, filter_mode='sos',dtype=dtype))
        else:
            reclist.append(subrec)
    return ChannelsAggregationRecording(reclist,recording.get_channel_ids())

def compute_notches_structs(recording,segment=0,span_smooth = 21, db_thr = 5,min_freq = 300, max_freq = 3000, spectum_resolution = 0.5, k_periodograms = 200):
    sr = recording.get_sampling_frequency()
    lts = recording.get_num_frames(segment)

    N = 2**np.ceil(np.log2(sr/spectum_resolution))
    samples_spectrum = int(min(k_periodograms *N ,lts))
    notches_dict = {}
    chIDs = recording.get_channel_ids()
    selected_fs = None
    for chi in range(recording.get_num_channels()):
        fs, pxx = welch(recording.get_traces(segment_index=segment,start_frame=0,end_frame=samples_spectrum, channel_ids=[chIDs[chi]])[:,0], 
            fs=sr,noverlap=0,window=windows.barthann(N), return_onesided=True)
        
        if selected_fs is None :
            selected_fs = np.logical_and(fs>=min_freq,fs<=max_freq)
        fs = fs[selected_fs]
        pxx = pxx[selected_fs]
        #pxx is a power spectral density 
        pxx_db = 10*np.log10(pxx)
        pxx_thr_db = medfilt(pxx_db,kernel_size = span_smooth) + db_thr


        low_freqs = int(np.maximum(fs.max(axis=0),min_freq) )           
        high_freqs = int(np.maximum(fs.max(axis=0),max_freq))
        pxx_thr_db[0:low_freqs] = pxx_thr_db[low_freqs]
        pxx_thr_db[high_freqs:-1] = pxx_thr_db[high_freqs]


        used_notches = []
        abs_db = []
        diff_db = []

        supra_thr = np.nonzero(pxx_db > pxx_thr_db)[0] 
        if len(supra_thr)>0:
            max_amp4notch = max(pxx_db- pxx_thr_db)
            temp_supra = np.nonzero(np.diff(supra_thr)>1)[0]
            inds_to_explore = np.concatenate([supra_thr[[0]], supra_thr[temp_supra+1]])

            if len(temp_supra)==0:
                sample_above = len(supra_thr)
            else:
                sample_above = np.concatenate([temp_supra[[0]], np.diff(temp_supra), len(supra_thr)-np.nonzero(supra_thr==inds_to_explore[-1])[0]+1])
            notch_idx = np.nonzero(sample_above[0:len(inds_to_explore)+1]>1)[0] #only uses sample_above(jj)>1
            used_notches = np.zeros([len(notch_idx),1])
            bw_notches = np.zeros([len(notch_idx),1])
            abs_db = np.zeros([len(notch_idx),1])
            diff_db = np.zeros([len(notch_idx),1])
            
            for i in range(len(notch_idx)):
                jj = notch_idx[i]
                iaux = np.argmax(pxx_db[inds_to_explore[jj]:inds_to_explore[jj]+sample_above[jj]-1])    #introduces alignment
                centre_sample = inds_to_explore[jj] + (sample_above[jj]-1)/2
                ind_max = iaux + inds_to_explore[jj] - 1
                if np.mod(centre_sample,1)==0.5 and (pxx_db[int(np.floor(centre_sample))]> pxx_db[int(np.ceil(centre_sample))]):
                    centre_sample = int(np.floor(centre_sample))
                else:
                    centre_sample = int(np.ceil(centre_sample))
                amp4notch = pxx_db[ind_max]-pxx_thr_db[ind_max]
                used_notches[i] = fs[centre_sample]
                bw_notches[i] = (fs[1]-fs[0])*sample_above[jj]*2*amp4notch/max_amp4notch
                abs_db[i] = pxx_db[ind_max]
                diff_db[i] = amp4notch

        bw_notches = np.zeros([len(used_notches),1])
        notches = {'Z': [],
        'P': [],
        'Ks': np.zeros([len(used_notches),1]),
        'freq': used_notches,
        'abs_db': abs_db,
        'diff_db': diff_db
        }

        for i in range(len(used_notches)):
            nfi = used_notches[i]
            w = nfi/(sr/2)
            if nfi<295:
                max_bw = 3
            else:
                max_bw = 5

            bw_notches[i] = min(max(bw_notches[i],1), max_bw)
            bw = bw_notches[i]/(sr/2)

            [b_notch,a_notch] = iirnotch(w, w/bw)
            [zi,pi,ki] = tf2zpk(b_notch,a_notch)
            notches['Z'].extend(zi)
            notches['P'].extend(pi)
            notches['Ks'][i] = ki
        notches['K'] = np.prod(notches['Ks'])
        notches['sos'] = zpk2sos(notches['Z'], notches['P'],notches['K'])
        notches['g']= 1
        notches['bw_notches'] = bw_notches
        notches_dict[chIDs[chi]] = notches
    return notches_dict
