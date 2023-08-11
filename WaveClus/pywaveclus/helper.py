import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from spikeinterface import ChannelSliceRecording
from scipy.signal import welch,windows

import numpy as np


def calculate_power_spectrum(recording, segment=0, spectum_resolution=0.5, k_periodograms=200):
    """
    Calculate the power spectrum of a recording. The power spectrum is calculated using Welch's method.
    We calculate the power_spectrum for the last segment of the recording.
    
    Args:
        recording: Recording object.
        segment: Index of the segment to calculate the power spectrum for (default: 0).
        spectum_resolution: Resolution of the power spectrum in Hz (default: 0.5).
        k_periodograms: Number of periodograms to average (default: 200).
        
    Returns:
        frequencies: Array of frequencies (in Hz).
        power_spectrum: Power spectrum.
            
    """    
    sampling_rate = recording.get_sampling_frequency()
    
    
    num_frames = recording.get_num_frames(segment)
    num_dft_points = int(2 ** np.ceil(np.log2(sampling_rate / spectum_resolution)))
    num_spectrum_samples = int(min(k_periodograms * num_dft_points, num_frames))
    channel_id = recording.get_channel_ids()[0]
    
    frequencies, power_spectrum = welch(
        recording.get_traces(
            segment_index=segment,
            start_frame=num_frames - num_spectrum_samples,
            end_frame=num_frames,
            channel_ids=[channel_id])[:, 0],
        fs = sampling_rate,
        noverlap=0,
        window=windows.barthann(num_dft_points),
        return_onesided=True
    )

    return frequencies, power_spectrum


def get_traces_for_time_range(recording, start_time, end_time, segment_index=0,ch_idx=0):
    """
    Retrieve traces for a specific time range from a recording.

    Args:
        recording: Recording object.
        start_time: Start time of the desired time range in seconds.
        end_time: End time of the desired time range in seconds.
        segment_index: Index of the segment to retrieve traces from (default: 0).

    Returns:
        time_points: Array of time points (in seconds) within the specified time range.
        traces: Traces within the specified time range.

    """
    
    fs = recording.get_sampling_frequency()  # Sampling frequency in Hz

    # Convert start and end time to frames
    start_frame = int(start_time * fs)
    end_frame = int(end_time * fs)

    # Retrieve traces for the specified time range
    traces = recording.get_traces(
        segment_index=segment_index, 
    start_frame=start_frame, 
    end_frame=end_frame,
    channel_ids=[ch_idx])[:, 0]

    # Create an array of time points (in seconds)
    time_points = np.arange(start_frame, end_frame) / fs

    return time_points, traces



