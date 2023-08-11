# Processing Pipeline

The following is a processing pipeline for analyzing neural recordings using MATLAB and Python. Here are the requirements for each:

### MATLAB Requirements  
- MATLAB version 9.6 (R2019a) or later
- Signal Processing Toolbox
- DSP System Toolbox
- Parallel Computing Toolbox
- MATLAB Parallel Server
- Polyspace Bug Finder
- [NPMK](https://github.com/BlackrockNeurotech/NPMK) for Blackrock recordings
- Neuroshare for Ripple recordings

### Python Requirements
You can install the required Python packages using the following command:
```
pip install -r requirements.txt
```

## Steps

1. **Parsing the Recording**
   - For Ripple recordings, use the `parse_ripple.m` script.
   - For Blackrock recordings, use the `parse_NSx.m` script.
   - Note: The Neuroshare library is required for Ripple recordings.

2. **Analyzing Power Spectrum and Calculating Notches**
   - Use the `new_check_lfp_power_NSX.m` script to generate figures for analyzing the power spectrum and calculating notches.
   - If you prefer to use Python for this step, there is a function implemented in Python called `filter_freq_peaks` that performs notch filtering. If you used the `new_check_lfp_power_NSX.m` script, you can set `load_mat_notches=False` in the Python pipeline.

3. **Configuration Setup**
   - Before starting the pipeline, make sure to set up the `config.Yaml` file in the `config` folder.
   - Update the paths in the configuration file to match the paths on your computer.

4. **Execution**
   - Open the `main.ipynb` notebook and follow the cells to execute the processing pipeline.
   - The notebook contains comments specifically for new Python users and provides details on how to use the `spikeinterface` library.


**Notes**
We need to change the bandpass filter from butterworth to epilettic


