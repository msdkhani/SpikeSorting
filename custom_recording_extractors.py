from typing import List, Union
from spikeinterface import BaseRecording, BaseRecordingSegment
from  spikeinterface.core.core_tools import read_binary_recording
import numpy as np
from scipy.io import loadmat
import warnings
from pathlib import Path
import yaml
import re

class NSXRecordingExtractor(BaseRecording):
    """
    RecordingExtractor for the NC5 format
    Parameters
    ----------
    folder_path: str or Path
        Path to the folder with the binary files
    channel_ids: list (optional)
        A list of channel ids
    Returns
    -------
    recording: NSXRecordingExtractor
        The recording Extractor
    """
    extractor_name = 'NSXRecordingExtractor'
    has_default_locations = False
    installed = True  # check at class level if installed or not
    is_writable = False
    mode = 'file'
    installation_mesg = ""  # error message when not installed
    
    def apply_gain(self, rec_segment, gain,dc):
        for i, timeseries in enumerate(rec_segment._timeseries):
            rec_segment._timeseries[i] = timeseries * gain + dc
        return rec_segment
    


    def get_channels_by_bundle(self):
        """
        Get a dictionary where each bundle contains a list of channels for the entire recording.
        Returns
        -------
        dict
            Dictionary where each bundle name maps to a list of dictionaries containing channel ID and label.
        """
        bundles = self.get_property("bundle")
        channel_ids = self.get_channel_ids()
        labels = self.get_property("label")
        channels_by_bundle = {}
        for channel_index, bundle_name in enumerate(bundles):
            channel_id = channel_ids[channel_index]
            label = labels[channel_index]
            channel_info = {"channel_id": channel_id, "label": label}
            if bundle_name in channels_by_bundle:
                channels_by_bundle[bundle_name].append(channel_info)
            else:
                channels_by_bundle[bundle_name] = [channel_info]
        return channels_by_bundle

    def __init__(self, channels, folder_path='.', probes_from = 'labels'):

        self._folder_path = Path(folder_path)
        metadata = loadmat(str(self._folder_path / 'NSx.mat'),variable_names='NSx', squeeze_me=True, simplify_cells=True)

        labels = []
        electrode_ID = []
        sr = None
        file_paths = []
        self._lts = np.inf
        self._ch2pos = {}
        det_channels = [] #channesl found
        probes = []
        gain = None
        dc = None
        bundles = []
        
        if channels is None:
            channels = [x['chan_ID'] for x in metadata['NSx'] if x['sr']==30000 and x['is_micro']==1]
        
        for ch in channels:
            #fist check the nc6
            info = list(filter(lambda x: (x['chan_ID'] == ch) and (x['ext'] != '.NC6'), metadata['NSx']))
            if len(info) == 0:
                info = list(filter(lambda x: (x['chan_ID'] == ch) and (x['ext'] == '.NC6'), metadata['NSx']))
            if len(info) == 0:
                assert info is not None, 'channel: {} not found'.format(ch)
            if len(info) == 0:
                warnings.warn('channel {} not found'.format(ch))
                continue
            det_channels.append(ch)
            info = info[0]
            if sr is None:
                sr = info['sr']
                gain = info['conversion']
                dc = info['dc']
            else:
                assert (sr == info['sr']), 'channels with different sampling rate'
                assert (gain == info['conversion']), 'channels with different sampling gain'
                assert(dc == info['dc']), 'channels with different dc offset'
            self._lts = np.min([self._lts, info['lts']]).astype(int)
            labels.append(info['label'])
            bundles.append(info['bundle'])
            electrode_ID.append(info['electrode_ID'])
            if probes_from == 'labels':
                probes.append(info['label'].split(' raw')[0])
            elif probes_from == 'channel/8':
                probes.append(np.floor((info['electrode_ID']- 1) / 8).astype(int))
            elif probes_from == 'channel/9':
                probes.append(np.floor((info['electrode_ID']- 1) / 9).astype(int)) #only for ripple systemswith the reference in each probe
            else:
                assert False, 'invalid probes_from input: {}'.format(labels)
            file = '{}{}'.format(info['output_name'], info['ext'])
            file_paths.append(str(self._folder_path / file))
            #self._timeseries.append(np.memmap(str(self._folder_path / file), dtype=np.int16,mode='r'))
            #self._ch2pos[ch] = len(self._timeseries)-1

        BaseRecording.__init__(self, sr, det_channels, np.int16)
        
        self.set_channel_gains(gain)
        self.set_channel_offsets(dc)
        self.annotate(is_filtered=False)
        self.set_property("electrode_ID", electrode_ID)
        self.set_property("bundle", bundles)
        self.set_property("label", labels)
        self.set_property("probe", probes)
        self.set_property("channel", det_channels)

        labels_wm =  [re.sub(r'\d*-\d*','',x) for x in labels]

        macros = {a:[i,b,0] for i,[a,b] in enumerate(zip(*np.unique(labels_wm,return_counts=True)))}
        nmacros = len(macros)
        positions = np.zeros((len(det_channels), 2))

        #just positions to plot
        r_plmicros = 8 #just example
        r_plmacros = r_plmicros * 5

        for i, m in enumerate(labels_wm):
            centerx = r_plmacros*np.cos(macros[m][0]*2*np.pi/nmacros)
            centery = r_plmacros*np.sin(macros[m][0]*2*np.pi/nmacros)
            positions[i,0] = centerx + r_plmicros*np.cos(macros[m][2]*2*np.pi/ macros[m][1])
            positions[i,1] = centery + r_plmicros*np.sin(macros[m][2]*2*np.pi/ macros[m][1])
            macros[m][2] = macros[m][2] + 1

        self.set_channel_locations(positions)
        rec_segment = MultiFileBinaryRecordingSegment(sr,file_paths, dtype=np.int16, file_offset=0)
        # multibly by conversion factor
        rec_segment = self.apply_gain(rec_segment, gain,dc)
        self.add_recording_segment(rec_segment)

        self._kwargs = {'folder_path': folder_path,
                        'channels': det_channels
                        }

 
class MultiFileBinaryRecordingSegment(BaseRecordingSegment):
    def __init__(self, sr,datfiles, dtype, file_offset):
        BaseRecordingSegment.__init__(self,sampling_frequency=sr)
        self._timeseries = []
        for file in datfiles:
            self._timeseries.append(read_binary_recording(file, num_chan=1, dtype=dtype, time_axis=0, offset=file_offset))
    def get_num_samples(self) -> int:
        """Returns the number of samples in this signal block
        Returns:
            SampleIndex: Number of samples in the signal block
        """
        return self._timeseries[0].shape[0]

    def get_traces(self,
                   start_frame: Union[int, None] = None,
                   end_frame: Union[int, None] = None,
                   channel_indices: Union[List, None] = None,
                   ) -> np.ndarray:
        
        if isinstance(channel_indices,slice):
            timeseries = self._timeseries[channel_indices]
        elif channel_indices is None:
            timeseries = self._timeseries
        else:
            timeseries = [self._timeseries[int(ix)] for ix in channel_indices.flatten()]
        return np.hstack([ch[start_frame:end_frame,:] for ch in timeseries])


with open("../config/config.yaml", 'r') as stream:
    config = yaml.safe_load(stream)

def load_notches_structs(recording):
    info = loadmat(config['mat_notches']['path'],
                   variable_names='process_info', squeeze_me=True, simplify_cells=True)['process_info']
    ids = [d['chID'] for d in info]
    sos = {}
    g = {}
    freq = {}
    chIDs = recording.get_channel_ids()
    notches_dict = {}
    for ch_id in chIDs:
        if ch_id in ids:
            sos = info[ids.index(ch_id)]['SOS'].copy(order='C')
            g = info[ids.index(ch_id)]['G']
            freq = info[ids.index(ch_id)]['freqs_notch']
        else:
            sos = None
            g = None
            freq = None
            
        notches_dict[ch_id] ={'sos': sos, 'g': g, 'freq': freq}
    return notches_dict

__version__ = '2.00'