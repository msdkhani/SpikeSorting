import numpy as np
import matplotlib.pyplot as plt
from helper import get_traces_for_time_range
from matplotlib.ticker import MultipleLocator
from spikeinterface import ChannelSliceRecording
from scipy.signal import medfilt
import os
import yaml
    
# Load the config object as a global variable
with open("../config/config.yaml", 'r') as stream:
    config = yaml.safe_load(stream)

# Define a global function to access the config object
def get_config():
    return config



#plot power spectrum in notches with raw data and bandpass filter and notches
def plot_notches(recording,recording_bp, recording_notched, pnotchinfo, fs,pxx,fs_raw,pxx_raw, min_freq=0, max_freq=10000):
    # Get channel ID
    channel_id = recording.get_channel_ids()[0]
    

    # Select frequencies within the desired range
    selected_frequencies = (fs >= 300) & (fs <= 3000)

    # Compute the thresholded power spectrum in decibels
    power_spectrum_db = 10 * np.log10(pxx)
    power_spectrum_db_raw = 10 * np.log10(pxx_raw)

    # Compute threshold power spectrum
    threshold_power_spectrum_db = medfilt(power_spectrum_db, kernel_size=21) + 5
    before_threshold_indices = np.where(selected_frequencies)[0].min()
    after_threshold_indices = np.where(selected_frequencies)[0].max()
    threshold_power_spectrum_db[:before_threshold_indices] = threshold_power_spectrum_db[before_threshold_indices]
    threshold_power_spectrum_db[after_threshold_indices:] = threshold_power_spectrum_db[after_threshold_indices]
    num_frames = recording.get_num_frames()
    sampling_rate = recording.get_sampling_frequency()
    
    end_time = num_frames/sampling_rate
    start_time = end_time - 60
    

    # Get timepoints and traces for the selected time range
    timepoints_raw, traces_raw = get_traces_for_time_range(recording, start_time, end_time, ch_idx=channel_id)
    timepoints_bandpass, traces_bandpass = get_traces_for_time_range(recording_bp, start_time, end_time, ch_idx=channel_id)
    timepoints_notched, traces_notched = get_traces_for_time_range(recording_notched, start_time, end_time, ch_idx=channel_id)

    # Select frequencies within the desired range for raw data
    selected_frequencies_raw = (fs_raw >= min_freq) & (fs_raw <= max_freq)

    # Create subplots
    fig, axs = plt.subplots(3, 1, sharex=False, figsize=(20, 15))

    # Plot power spectra
    axs[0].plot(fs_raw[selected_frequencies_raw], power_spectrum_db_raw[selected_frequencies_raw], label='Power Spectrum raw data', color='royalblue',linewidth=1)
    axs[0].plot(fs[selected_frequencies_raw], power_spectrum_db[selected_frequencies_raw], label='Power Spectrum BandPass Filter', color='red',linewidth=1)
    axs[0].plot(fs[selected_frequencies_raw], threshold_power_spectrum_db[selected_frequencies_raw], label='Threshold', color='orchid',linewidth =2)
    #if pnotchinfo[channel_id]:
    #   axs[0].scatter(pnotchinfo[channel_id]['freq'], pnotchinfo[channel_id]['abs_db'], color='green', label='Notches')
    if pnotchinfo[channel_id]:
        for notch in pnotchinfo[channel_id]['freq']:
            axs[0].axvline(x=notch, color='r', linestyle='--',linewidth=0.7)
    # Set y-axis limits for power spectra
    axs[0].set_ylim([-50, 50])


    # Plot raw traces
    axs[1].plot(timepoints_raw, traces_raw, color='royalblue', linewidth=0.3, label=f'Channel {channel_id}')

    # Plot bandpass-filtered traces
    axs[2].plot(timepoints_bandpass, traces_bandpass, color='royalblue', linewidth=0.3, label=f'Channel {channel_id}')
    axs[2].plot(timepoints_notched, traces_notched, color='red', linewidth=0.3, label=f'Channel {channel_id}')

    # Set x-axis limits
    axs[0].set_xlim(0, 10000)
    axs[1].set_xlim(start_time, end_time)
    axs[2].set_xlim(start_time, end_time)

    # Set x and y labels
    axs[0].set_xlabel('Frequency (Hz)')
    axs[0].set_ylabel('Power Spectrum (dB)')
    axs[1].set_xlabel('Time (s)')
    axs[1].set_ylabel('Amplitude')
    axs[2].set_xlabel('Time (s)')
    axs[2].set_ylabel('Amplitude')

    # Add legends
    axs[0].legend(loc='upper right')
    axs[1].legend(loc='upper right')
    n_notched= len(pnotchinfo[channel_id]['freq'])
    fig.suptitle(f'Spectrum of Channel: {channel_id} Session check. Sr: {sampling_rate}, #Notches:{n_notched}, Filtered: Band-Psss: 300-3000', fontsize=14, fontweight='bold',y=0.9)
    
    config = get_config()
    
    PATH = config['paths']
    OUTPUT = PATH['OUTPUT_FOLDER']

    # Save the figure
    output_folder = config['paths']['OUTPUT_FOLDER']
    project_name = config['project']['name']
    save_dir = os.path.join(output_folder, project_name)
    save_dir = os.path.join(save_dir, 'plots')
    save_dir = os.path.join(save_dir, 'Spectrum')

    os.makedirs(save_dir, exist_ok=True)
    save_path = save_dir+ '/'+f'spectrum_check_CH{channel_id}.png'
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)  # Close the figure to free up memory

    print(f"Figure saved: {save_path}")
    
    # Display the plot
    #plt.show()
    
    
    

def plot_bundle_traces(recording_notched, notched_info, start_time=10, end_time=130, bundle_dict=[]):
    """
    Function to plot traces from a recording based on channel bundles.

    Args:
        recording_bp: Bandpass filtered recording object.
        recording_notched: Notched recording object.
        notched_info: Dictionary containing notch filter information for each channel.
        start_time: Start time of the plot range in seconds.
        end_time: End time of the plot range in seconds.
        bundle_dict: Dictionary containing channel bundles, where each bundle is a list of channel IDs.
    """
    
     # Retrieve custom settings from the global config object
    config = get_config()
    
    PATH = config['paths']
    OUTPUT = PATH['OUTPUT_FOLDER']
    
    start_time = config['plot_continuous']['start_time']
    end_time = config['plot_continuous']['end_time']
    factor_thr = config['plot_continuous']['factor_thr']
    w_pre = config['plot_continuous']['w_pre']
    w_post = config['plot_continuous']['w_post']
    min_ref_per = config['plot_continuous']['min_ref_per']
    
    # Customize plot appearance
    plt.style.use('seaborn-whitegrid')  # Apply a modern style

    # Iterate over bundles
    for bundle_name in bundle_dict:
        bundle = bundle_dict[bundle_name]
        
        # Create the figure and subplots for the current bundle
        fig, axs = plt.subplots(len(bundle), 1, sharex=True, figsize=(15, 8))

        # Check if there is only one channel in the bundle
        if len(bundle) == 1:
            axs = [axs]  # Convert single AxesSubplot to a list

        # Iterate over channels in the bundle
        for idx, channel_id in enumerate(bundle):
            # Create a sub-recording for the current channel
            #sub_recording_bandpass = ChannelSliceRecording(recording_bp, [channel_id])
            label = channel_id['label']
            channel_id = channel_id['channel_id']
            
            sub_recording_notched = ChannelSliceRecording(recording_notched, [channel_id])
            

            timepoints, traces= get_traces_for_time_range(sub_recording_notched, start_time, end_time, ch_idx=channel_id)
            sr = sub_recording_notched.get_sampling_frequency()
            

            ref = np.floor(min_ref_per * sr / 1000)
            thr = factor_thr * np.median(np.abs(traces)) / 0.6745
            thrmax = 10 * thr
            
            xaux = np.where((traces[w_pre + 1: -w_post - 1] < -thr) & (np.abs(traces[w_pre + 1: -w_post - 1]) < thrmax))[0] + w_pre


            num_spikes = np.count_nonzero(np.diff(xaux) > ref) # Number of spikes detected
            # Plot the downsampled data on the corresponding subplot
            axs[idx].plot(timepoints, traces, color='royalblue', linewidth=0.3, label=f'Channel {channel_id}')
            

            axs[idx].axhline(-thr, color='red', linewidth=0.6)

            axs[idx].set_ylabel(f'Channel {channel_id}', fontweight='bold', fontsize=10)
            axs[idx].set_title(f'{label}, Notches Applied : {len(notched_info[channel_id]["freq"])} , # Spikes : {num_spikes}, Threshold: {thr}', fontweight='bold', fontsize=10)

            # Customize y-axis ticks and limits
            axs[idx].set_ylim([-100, 50])
            axs[idx].tick_params(axis='both', which='both', labelsize=8)
            # Set X-axis tick locator
            axs[idx].xaxis.set_major_locator(MultipleLocator(10))

        # Set the x-label for the last subplot
        axs[-1].set_xlabel('Time (s)', fontweight='bold', fontsize=12)
        # Set the X-axis limits
        for ax in axs:
            ax.set_xlim(start_time, end_time)

        # Adjust the spacing between subplots
        plt.tight_layout()

        # Add legend outside the plot
        #fig.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')

        # Set the title of the figure as the bundle name
        fig.suptitle(f'Bundle {bundle_name}', fontsize=14, fontweight='bold',y=1.02)



        # Save the figure
        output_folder = config['paths']['OUTPUT_FOLDER']
        project_name = config['project']['name']
        save_dir = os.path.join(output_folder, project_name)
        save_dir = os.path.join(save_dir, 'plots')
        save_dir = os.path.join(save_dir, 'bundles')
        

        os.makedirs(save_dir, exist_ok=True)
        save_path = save_dir+ '/'+ f'bundle_{bundle_name}.png'
        plt.savefig(save_path, bbox_inches='tight')
        plt.close(fig)  # Close the figure to free up memory

        print(f"Figure saved: {save_path}")
        


