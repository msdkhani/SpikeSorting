import numpy as np
import itertools
from spikeinterface.comparison import compare_two_sorters
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from spikeinterface.postprocessing  import WaveformPrincipalComponent, compute_spike_amplitudes
from spikeinterface.qualitymetrics import compute_quality_metrics, calculate_pc_metrics
import spikeinterface as si
from os import listdir
from joblib import Parallel, delayed, cpu_count

sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})

leicolors_list = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.620690, 0.0, 0.0], [0.413793, 0.0, 0.758621],
                  [0.965517, 0.517241, 0.034483], [0.448276, 0.379310, 0.241379], [1.0, 0.103448, 0.724138],
                  [0.545, 0.545, 0.545], [0.586207, 0.827586, 0.310345], [0.965517, 0.620690, 0.862069],
                  [0.620690, 0.758621, 1.]] #silly name, just colors that look different enough
leicolors = lambda x: leicolors_list[x % len(leicolors_list)]
from scipy.stats import gaussian_kde
import matplotlib 

def split_sorting_by(sorting, property, outputs='dict',rec=None):
        assert outputs in ('list', 'dict')
        values = sorting.get_property(property)
        if values is None:
                
            if (rec is not None 
                and 'channel' in rec.get_property_keys()
                and  'channel' in sorting.get_property_keys()):
                values = rec.get_property('probe')
                rchannels = rec.get_property('channel')
                schannels = sorting.get_property('channel')
                values = np.array([values[np.nonzero(rchannels==int(s))[0]][0] for s in schannels])
            else:
                raise ValueError(f'property {property} is not set')
            
        if outputs == 'list':
            sortings = []
        elif outputs == 'dict':
            sortings = {}
        
        for value in np.unique(values):
            inds, = np.nonzero(values == value)
            sub_sort = sorting.select_units(sorting.unit_ids[inds], 
                renamed_unit_ids=np.arange(1, 1+len(inds)))
            if outputs == 'list':
                sortings.append(sub_sort)
            elif outputs == 'dict':
                if isinstance(value, str) and value.isnumeric():
                    value = int(value)
                sortings[value] = sub_sort
        return sortings



def same_probe_collisions(sorting_all,recording, prop_thr=-1, delta_time=0.4, outputfolder = '.'):
    
    
    results = {}
    sortings_probe = split_sorting_by(sorting_all, 'probe', outputs='dict',rec=recording)
    # here parallel can be added (FC)
    for probe,sorting_probe in sortings_probe.items():
        result = []
        sortings = split_sorting_by(sorting_probe,'channel',outputs = 'dict')
        pairs = itertools.combinations(sortings.keys(),2)
        for ch1, ch2 in pairs:
            comp = compare_two_sorters(sortings[ch1], sortings[ch2], delta_time=delta_time)
            for label2, col in comp.match_event_count.items():
                for label1, match_count in col.iteritems():
                    coefb = match_count/comp.event_counts2[label2]
                    coefa = match_count/comp.event_counts1[label1]
                    if  coefb >= prop_thr or coefa >= prop_thr:
                        #if  coefa >= prop_thr:
                        result.append({'Cl-Ch_A':'{}-{}'.format(label2,ch2),'Cl-Ch_B':'{}-{}'.format(label1,ch1),'coef':coefa})
                        #if  coefb >= prop_thr:
                        result.append({'Cl-Ch_A':'{}-{}'.format(label1,ch1),'Cl-Ch_B':'{}-{}'.format(label2,ch2),'coef':coefb})

        if result:
            data = pd.pivot_table(pd.DataFrame(result), index=['Cl-Ch_A'], columns=['Cl-Ch_B'], fill_value=0)['coef']
            results[probe] = data.iloc[sort_Cl_Ch(data.index.values),sort_Cl_Ch(data.columns.values)]

    for i,(probe, matrix) in enumerate(results.items()):
        fig = plt.figure(dpi=100)
        ax = sns.heatmap(matrix, linewidths=0.5, cbar = True)
        ax.set_title(label='{}: $ \\| (A \\cap B) \\| / \\|B \\| $'.format(probe))
        plt.tight_layout()
        plt.savefig(outputfolder/Path('collision_{}.png'.format(probe)))
        print(f'Collisions derected in probe:{probe}')
        plt.close(fig)
    return results

def sort_Cl_Ch(series):
    new_order = []
    for x in series:
        parts = x.split('-')
        new_order.append(int(parts[1])*2000+int(parts[0]))
    return np.argsort(new_order)


def plot_sorting_results(output_folder, channels=None, **kwargs):
    '''
    For each channel in the dictionaty, it Saves figures with mean spikes waveform comparisons, spikes waveforms for each unit, ISI histograms and time activity
    of each unit.

    Parameters
    ----------
    output_folder: str
        Folder where the figures will be saved.
    channels: list or None (default)
        If is not None, channesl that will be used.
    kwargs: arguments of plot_channels_from_wf
    '''
    available_channels = [int(x[6:]) for x in listdir(output_folder) if x.startswith('wf_ch_')]
    if channels is None:
        channels = available_channels
    else:
        channels = set(channels).intersection(available_channels)

    delayed_funcs = [delayed(plot_channels_from_wf)(output_folder/Path(f'wf_ch_{c}'),output_folder,**kwargs) for c in channels]
    parallel_pool = Parallel(n_jobs = min(cpu_count(), len(channels)))

    parallel_pool(delayed_funcs)


def plot_channels_from_wf(f, output_folder, nspikes=3000, bin_step=1, nbins=100, time_pixes=200, prefix='units',dpi=100):
    '''
    Plot figure showing waveforms, ISI and number of spikes.
    Parameters
    ----------
    f: folder of a WaveformExtractor
    output_folder: str
        Folder where the figures will be saved.
    channels: list or None (default)
        If is not None, channesl that will be used.
    nspikes: int, default: 3000
    bin_step: int, default: 1
        ISI bin in ms.
    nbins: int, default: 100
    time_pixes: int, default: 200
        number of horizontal pixels on the activity plot
    prefix: str, default: 'units'
        Prefix of the created figures. Can help for saving after editing the sorting
    '''
    
    we=si.WaveformExtractor.load_from_folder(f)
    x4plot = np.arange(1, we.nsamples+1)
    sns.set_style('seaborn-whitegrid')
    plt.rcParams.update({'font.size': 18})        
    from matplotlib.colors import LogNorm

    sorting = we.sorting
    recording = we.recording

    time_grid = np.linspace(0, recording.get_num_frames() / (sorting.get_sampling_frequency() * 60), time_pixes)
    pix2min = recording.get_num_frames() / (sorting.get_sampling_frequency() * 60 * time_pixes)

    srk = sorting.get_sampling_frequency() / 1000
    units = sorting.get_unit_ids()
    spikes = {}
    means = {}
    stds = {}
    for i, u in enumerate(units):
        spikes[u] = we.get_waveforms(u)[:,:,0]
        means[u] = np.mean(spikes[u], 0)
        stds[u] = np.std(spikes[u], 0)

    time_den = []
    ulims = {}  # for axis
    if len(units) > 4:
        nfigures = 1 + int(np.ceil((len(units)-4)/5))
    else:
        nfigures=1
    for fig_num in range(nfigures):
        (fig, (axs1, axs2, axs3)) = plt.subplots(3, 5, figsize=(24, 12), dpi=dpi)
        if fig_num == 0:
            units2p = [units[u] for u in range(min(4, len(units)))]
        else:
            units2p = [units[u] for u in range(4 + 5 * (fig_num - 1), min(4 + 5 * fig_num, len(units)))]

        for a in axs1:
            a.autoscale(False)

        if fig_num == 0:
            axs1[0].grid('both')
            for i, u in enumerate(units):
                qs = np.quantile(spikes[u], [0.01, 0.99], axis=0)
                ulims[u] = [qs[0].min() * 1.1, qs[1].max() * 1.1]
                axs1[0].plot(x4plot,means[u], color=leicolors(i), linewidth=3)
                sptrain = sorting.get_unit_spike_train(u) / (sorting.get_sampling_frequency() * 60)
                if len(sptrain) > 1:
                    try:
                        time_den.append(gaussian_kde(sptrain).pdf(time_grid))
                    except np.linalg.LinAlgError:
                        aux = np.zeros_like(time_grid)
                        aux[np.argmin(np.abs(time_grid - sptrain[0]))] = 1
                        time_den.append(aux)
                else:
                    aux = np.zeros_like(time_grid)
                    aux[np.argmin(np.abs(time_grid - sptrain))] = 1
                    time_den.append(aux)
            axs1[0].set_title('Total spikes #{}'.format(sum([w.shape[0] for w in spikes.values()])))
            if len(units) > 0:
                we.load_extension('quality_metrics').get_data()['l_ratio'].plot.bar(ylabel = 'Lratio', xlabel = 'Cluster', ax=axs3[0],rot=0)
                axs2[0].imshow(np.vstack(time_den), cmap=plt.cm.inferno, aspect='auto',interpolation='none')
                axs2[0].hlines(np.arange(len(units)) + 0.5, 0, len(time_den[0]), color='k', linewidth=4)
                axs2[0].set_yticks(np.arange(len(units)))
                axs2[0].set_yticklabels(['Cl: {}'.format(u) for u in units], fontsize=10)
                axs2[0].set_xticks(axs2[0].get_xticks())
                axs2[0].set_xticklabels(['{:.1f}'.format(xt*pix2min) for xt in axs2[0].get_xticks()])
                axs2[0].set_xlim(0, len(time_grid))
                axs2[0].tick_params(axis='y', which='both', grid_linestyle='None')
                axs2[0].set_xlabel('Time (min)')
                axs2[0].set_ylabel('Presence Plot')
                ylims = [min([y[0] for y in ulims.values()]), max([y[1] for y in ulims.values()])]  # general
        for i, u in enumerate(units2p):
            ls = means[u].shape[0]
            if fig_num == 0:
                axi = i+1
                axs1[0].set_ylim(ylims)
                axs1[0].set_xlim([1, ls])
                axs1[0].set_xlabel('Sample')
            else:
                axi = i
            allunits_i = np.where(units==u)[0][0]
            index2plot = np.random.choice(np.arange(spikes[u].shape[0]), size=min(nspikes,spikes[u].shape[0]), replace=False) 
            axs1[axi].plot(x4plot, spikes[u][index2plot,:].T, color=leicolors(allunits_i), alpha=0.3)
            axs1[axi].plot(x4plot, means[u], color='k', linewidth=3)
            axs1[axi].plot(x4plot, means[u] + stds[u], color='k', linewidth=1, linestyle='--')
            axs1[axi].plot(x4plot, means[u] - stds[u], color='k', linewidth=1, linestyle='--')
            axs1[axi].grid('both')
            axs1[axi].set_title('Cluster: {} #{}'.format(u, spikes[u].shape[0]))

            hpixels = np.round(ls * 2 / 3).astype(int)
            grid = np.linspace(ulims[u][1], ulims[u][0], hpixels)
            ks = []
            for a in spikes[u].T:
                if len(a) > 1:
                    try:
                        ks.append(gaussian_kde(a).pdf(grid))
                    except np.linalg.LinAlgError:
                        aux = np.zeros_like(grid)
                        aux[np.argmin(np.abs(grid - a[0]))] = 1
                        ks.append(aux)
                else:
                    aux = np.zeros_like(grid)
                    aux[np.argmin(np.abs(grid - a))] = 1
                    ks.append(aux)
                    
            matden = np.vstack(ks).T
            #, vmax=np.quantile(matden.flatten(),0.98) # to "improve" heatmap
            axs2[axi].imshow(matden, cmap=plt.cm.inferno, aspect='auto',interpolation='none')
            yticks = np.linspace(0, hpixels-1, 5, dtype=int)
            xticks = np.arange(-1, ls, 20, dtype=int)
            axs2[axi].set_yticks(yticks)
            axs2[axi].set_yticklabels(['{:.1f}'.format(g) for g in grid[yticks]])
            axs2[axi].set_xticks(xticks)
            axs2[axi].set_xticklabels(['{}'.format(g+1) for g in xticks])

            axs2[axi].tick_params(axis=u'both', which=u'both', grid_linestyle=':')

            times_diff = np.diff(sorting.get_unit_spike_train(u) / srk)
            axs3[axi].autoscale(False, axis='x')
            axs3[axi].hist(times_diff, bins=np.arange(0, nbins, bin_step),rwidth=1,linewidth=0)
            axs3[axi].set_title('{} in < 3ms'.format(np.count_nonzero(times_diff < 3)))
            axs3[axi].set_xlabel('ISI (ms)')
            axs3[axi].set_xlim([0, bin_step * nbins])

            axs1[axi].set_ylim(ylims)
            axs1[axi].set_xlim([1, ls])
            axs1[axi].set_xlabel('Sample')

        axs1[0].set_ylabel('Amplitude (uV)')

        output_folder = Path(output_folder)
        figname = '{}_ch{}_{}'.format(prefix, recording.get_channel_ids()[0], fig_num)
        fig.suptitle(str(output_folder / figname))
        #fig.tight_layout()
        output_folder.mkdir(exist_ok=True)
        fig.savefig(str(output_folder / (figname+'.png')))
        plt.close(fig)


def create_waveform_extractors_by_channel(sorting_all, recording, output_folder, recompute_chs=[]):
    """
    IT creates waveforms extractos for each channel and compute metrics
    ----------
    sorting_all: SortingExtractor
        Channel should be a propierty.
    recordings_ch: RecordingExtractor
        Channel should be a propierty.
    output_folder: str
        output folder
    Returns
    -------
    recording: NSXRecordingExtractor
        The recording Extractor
    """

    # it splits the recording and sorting by channel , creating 2 dictionaries
    sortings_ch = split_sorting_by(sorting_all, 'channel', outputs='dict',rec=recording)
    recordings_ch = recording.split_by('channel', outputs='dict')
    recompute_chs = [f'{c}' for c in recompute_chs] #pass to string because the propierty could change its type after loading
    waveforms_ch = {}
    # here parallel can be added (FC)
    for channel, sorting in sortings_ch.items():
        if recompute_chs and not any([f'{channel}'==c for c in recompute_chs]):
            continue
        waveform_folder = output_folder/Path(f'wf_ch_{channel}')
        we = si.WaveformExtractor.create(recordings_ch[channel], sorting, waveform_folder,remove_if_exists=True)
        
        we.set_params(ms_before=0.65, ms_after=1.46, max_spikes_per_unit=9999999,return_scaled=True)
        we.run_extract_waveforms(n_jobs=4,chunk_memory='700M',verbose=False)

        pca = WaveformPrincipalComponent(we)
        pca.set_params(n_components=5, mode='by_channel_local')
        pca.run()
        compute_quality_metrics(we)
        calculate_pc_metrics(pca)
        compute_spike_amplitudes(we, peak_sign='neg')
        waveforms_ch[channel] = we
    return waveforms_ch

def load_waveforms_extractors(OUTPUT_FOLDER, channels=None):
    waveforms_ch = {}
    if channels is None:
        channels = [int(x[6:]) for x in listdir(OUTPUT_FOLDER) if x.startswith('wf_ch_')]
    for channel in channels:
        waveforms_ch[channel] = si.WaveformExtractor.load_from_folder(OUTPUT_FOLDER/Path(f'wf_ch_{channel}'))
    return waveforms_ch