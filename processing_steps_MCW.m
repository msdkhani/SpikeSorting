%% FIRST: Move to the folder where the data is stored.
function processing_steps_MCW(which_system_micro,fast_analysis,nowait,micros,remove_ref_chs,do_power_plot,notchfilter,do_sorting,do_loop_plot,extract_events,extra_stims_win)

clearvars
clc
% filenames={'sHiCikEPkEMkGm_20161018-160552-001'};
%  which_system_micro = 'BRK'; % 'BRK' or 'RIP'
if ~exist('which_system_micro','var')|| isempty(which_system_micro),  which_system_micro = 'RIP'; end % 'BRK' or 'RIP'
if ~exist('fast_analysis','var')|| isempty(fast_analysis), fast_analysis = 0; end % fast_analysis = 1;
if ~exist('nowait','var')|| isempty(nowait), nowait = true; end % nowait = false;
if ~exist('micros','var')|| isempty(micros), micros = true; end % micros = false;
if ~exist('remove_ref_chs','var')|| isempty(remove_ref_chs), remove_ref_chs = []; end % remove_ref_chs = [265,274,297,306,329,338];
if ~exist('do_power_plot','var')|| isempty(do_power_plot), do_power_plot = 1; end 
if ~exist('notchfilter','var')|| isempty(notchfilter), notchfilter=1; end 
if ~exist('do_sorting','var')|| isempty(do_sorting), do_sorting = 1; end 
if ~exist('do_loop_plot','var')|| isempty(do_loop_plot), do_loop_plot = 1; end 
if ~exist('extract_events','var')|| isempty(extract_events), extract_events = 1; end 
if ~exist('extra_stims_win','var')|| isempty(extra_stims_win), extra_stims_win = 0; end 

RIP_hours_offset = 5;
NSPs={'_NSP-1';'_NSP-2'};
USE_PHOTODIODE = 1;
ch_photo_BRK = 1257;
ch_photo_RIP = 10241;
step_pic = 15;


%%
% Define the path where the codes emu repository is located
[~,name] = system('hostname');
if contains(name,'BEH-REYLAB'), dir_base = '/home/user/share/codes_emu'; 
elseif contains(name,'TOWER-REYLAB') || contains(name,'RACK-REYLAB'),
%     current_user = 'sofiad';  % replace with appropriate user name  
    current_user = getenv('USER');    dir_base = sprintf('/home/%s/Documents/GitHub/codes_emu',current_user); 
elseif contains(name,'NSRG-HUB-15446'), dir_base = 'D:\codes_emu'; % Hernan's desktop
end

addpath(dir_base);
custompath = reylab_custompath({'wave_clus_reylab','NPMK','codes_for_analysis','mex','useful_functions','neuroshare' });


if ~exist('dirmain','var')
    dirmain = pwd;
end

% addpath(genpath([dirmain '_pic'])) % srtimuli folder
set(groot,'defaultaxesfontsmoothing','off')
set(groot,'defaultfiguregraphicssmoothing','off')
set(groot,'defaultaxestitlefontsizemultiplier',1.1)
set(groot,'defaultaxestitlefontweight','normal')
%%
% exp_type = 'RSVPSCR';  phase=[];
exp_type = 'RSVP_online';  phase=[];
addpath(genpath([dirmain filesep 'pics_used']))
cd(dirmain)
if micros
    ftype = 'ns5';
else
    if strcmp(which_system_micro,'RIP')
        ftype = 'nf3';
    elseif strcmp(which_system_micro,'BRK')
        ftype = 'ns3';
    end
end

if strcmp(which_system_micro,'BRK')
    fname_prefix = '*.';
elseif strcmp(which_system_micro,'RIP')
    fname_prefix = '*_RIP.';
end

if ~exist('filenames','var')
    A=dir([fname_prefix ftype]);
    if isempty(A)
        error('There are no %s files in this folder',ftype);
    else
        filenames = {A.name};
    end
else
    fprintf('variable filenames already exists and is equal to %s\n',filenames{:})
end
%%
max_memo_GB = [];

if strcmp(which_system_micro,'BRK')
    if fast_analysis
        if ~exist('bundle_labels','var')
            parse_NSx(filenames,[1 2],max_memo_GB);
        else
            parse_NSx(filenames,[1 2],max_memo_GB,bundle_labels);
        end
    else
        if all(cellfun(@(x) ~contains(x,NSPs{1}) && ~contains(x,NSPs{2}),filenames))
            parse_NSx(filenames,[],max_memo_GB);
            filename_NEV = [filenames{1}(1:end-3) 'nev'];
        else
            parse_NSx(filenames(cellfun(@(x) contains(x,NSPs{2}),filenames)),2,max_memo_GB);
            ind_file=cellfun(@(x) contains(x,NSPs{2}),filenames);
            filename_NEV = [filenames{ind_file(1)}(1:end-3) 'nev'];
        end
    end
elseif strcmp(which_system_micro,'RIP')
    %         parse_ripple(filenames,remove_ref_chs,bundle_labels)
    parse_ripple(filenames,RIP_hours_offset,remove_ref_chs)
    filename_NEV = [filenames{1}(1:end-3) 'nev'];
end
clear filenames

%%
if ~exist('channels','var')
    channels=[];
    load('NSx','NSx');
    if micros
        AA = {NSx(arrayfun(@(x) (strcmp(x.unit,'uV') && x.sr==30000),NSx)).chan_ID};
    else
        AA = {NSx(arrayfun(@(x) (x.sr==2000),NSx)).chan_ID};
    end

    for i=1:length(AA)
        channels(i)=double(AA{i});
    end
end
% channels = [1:246];
%% opens parallel pool
poolobj = gcp('nocreate');
if isempty(poolobj) % If already a pool, do not create new one.
    parpool
end

%% power spectrum and potential notches
if do_power_plot
%     new_check_lfp_power_NSX(channels,'parallel',true)
    new_check_lfp_power_NSX(channels,'parallel',true,'with_NoNotch',true)
    disp('power spectra DONE')
end
%% plot continuous data (spike filtered) for each macro (we dont really gain time by running it in parallel, at least when plotting just 2 mins per channel)

% sequential
filt_order=4;
tic
Plot_continuous_bundles(channels,notchfilter,filt_order,'neg')
toc
disp('plot continuous DONE')
%% read pulses STEP 1: open NEV file (to read pulses sent with the DAQ)
if extract_events && ~fast_analysis && micros
    %ONLY FOR SINGLE FILE AT THE MOMENT
    if strcmp(which_system_micro,'BRK')
        openNEV(['./' filename_NEV],'report','noread','8bits')
    else
        % CHECK processing NEV
    end
    
    if strcmp(exp_type,'RSVP_online')
        if USE_PHOTODIODE
            extract_events_rsvpscr_BCM_photoonly_online3(which_system_micro,nowait)
        else
            if strcmp(which_system_micro,'RIP')
                extract_events_rsvpscr_ripple_EMU_online3(filename_NEV); % CHECK matlab times
            end
        end
    end
end
%% spike detection
if micros
    if fast_analysis || nowait
        neg_thr_channels = channels;
        pos_thr_channels = [];
    else
        ch_temp = input(sprintf('Currently, channels = %s.\nIf you want to keep it like that, press enter.\nOtherwise, enter the new vector and press enter  ',num2str(channels)));
        if ~isempty(ch_temp)
            channels = ch_temp;
        end

        neg_thr_channels = input('Enter the vector with neg_thr_channels and press enter. Press enter to use all channels ');
        if isempty(ch_temp)
            neg_thr_channels = channels;
            pos_thr_channels = [];
        else
            pos_thr_channels = input('Enter the vector with pos_thr_channels and press enter. Press enter for empty array ');
        end
    end

    parallel = true;
    % CHECK VALUES IN COMPARISON WITH SET_PARAMETERS.TXT
    param.detect_order = 4;
    param.sort_order = 2;
    param.detect_fmin = 300;
    param.sort_fmin = 300;
    param.stdmin = 5;
    param.stdmax = 50;                     % maximum threshold for detection
    param.ref_ms = 1.5;
    % param.preprocessing = false;
    param.preprocessing = true;

    tic
    param.detection = 'neg';
    % neg_thr_channels=[2097:2104];
    if ~isempty(neg_thr_channels)
        Get_spikes(neg_thr_channels,'parallel',parallel,'par',param);
    end

    % pos_thr_channels=[];
    param.detection = 'pos';
    if ~isempty(pos_thr_channels)
        Get_spikes(pos_thr_channels,'parallel',parallel,'par',param);
    end

    param.detection = 'both';
    % both_thr_channels=[65 67:74 76:80 82:86 88:96 113:122 126 128];
    both_thr_channels=setdiff(setdiff(channels,neg_thr_channels),pos_thr_channels);
    if ~isempty(both_thr_channels)
        Get_spikes(both_thr_channels,'parallel',parallel,'par',param);
    end
    toc
    disp('spike detection DONE')

    %% build data structure (grapes) and compute best responses for multiunits
    if do_loop_plot  && ~fast_analysis
        muonly = 'y';
    
        if strcmp(exp_type,'RSVPSCR')
            skip = 0; ons_ind =0;effect_rows_2=0;
            do_structure_mu_BCM_online3(channels,exp_type)
            rankfirst=1; ranklast=15;    
            loop_plot_best_responses_BCM_rank(channels,muonly,rankfirst,ranklast,effect_rows_2,1,3)
        end
    
    
        % same as before in case there are channels where more than the "best" 15 responses are needed
        channels_more = channels;
        rankfirst=16;
        for n_win = 1:extra_stims_win
            if strcmp(exp_type,'SCR') || strcmp(exp_type,'RSVPSCR')
                ranklast=rankfirst+step_pic-1;
                loop_plot_best_responses_BCM_rank(channels_more,muonly,rankfirst,ranklast,effect_rows_2,1,3)
            end
            rankfirst = rankfirst + step_pic;
        end
    %     disp('plot best responses DONE')
    end
    %% sorting
    if do_sorting
        param.min_clus = 15;
        param.max_spk = 30000;
        param.mintemp = 0.00;                  % minimum temperature for SPC
        param.maxtemp = 0.251;                 % maximum temperature for SPC
        param.tempstep = 0.01;
        
        Do_clustering(channels,'parallel',true,'make_times',true,'make_plots',false,'par',param)
    
        disp('spike sorting DONE')
    %%
        % build data structure (grapes) and compute best responses for single units
        if do_loop_plot  && ~fast_analysis
            clustered_channels = channels;
            muonly = 'n';
    
            if strcmp(exp_type,'SCR') || strcmp(exp_type,'RSVPSCR')
                rankfirst=1; ranklast=15;
            else
                rankfirst=1; ranklast=10;
            end
    
            if strcmp(exp_type,'RSVPSCR')
                do_structure_sorted_BCM_online3(clustered_channels)
    %             loop_plot_best_responses_BCM_rank(clustered_channels,muonly,rankfirst,ranklast,effect_rows_2,1,3)
            end
    
            channels_more_clus = clustered_channels;
            rankfirst=16;
            for n_win = 1:extra_stims_win
                if strcmp(exp_type,'SCR') || strcmp(exp_type,'RSVPSCR')
                    ranklast=rankfirst+step_pic-1;
    %                 loop_plot_best_responses_BCM_rank(channels_more_clus,muonly,rankfirst,ranklast,effect_rows_2,1,3)                
                end
                rankfirst = rankfirst + step_pic;
            end
    %         disp('plot best responses DONE')
        end
    end
    plot_grapes_as_online(1)
    disp('plot best responses DONE')
    Do_clustering(channels,'parallel',true,'make_times',false,'make_plots',true,'par',param)
end
%%
set(groot,'defaultaxesfontsmoothing','remove')
set(groot,'defaultfiguregraphicssmoothing','remove')
set(groot,'defaultaxestitlefontsizemultiplier','remove')
set(groot,'defaultaxestitlefontweight','remove')

rmpath([dirmain filesep 'pics_used'])
cellfun(@(x) rmpath(genpath(strjoin(x, filesep))),paths2add);
