function Plot_continuous_bundles(channels,notchfilter,filt_order,detect,tmin)
% Plots the continuous data and does a first detection of spikes for each
% channel. 
% detname = 1; 
if ~exist('filt_order','var')|| isempty(filt_order),     filt_order = 4; end
if ~exist('notchfilter','var')|| isempty(notchfilter),     notchfilter = 0; end
if ~exist('detect','var')|| isempty(detect),     detect = 'both'; end  %'pos','neg','both' detname = 0;
if ~exist('tmin','var')|| isempty(tmin),     tmin = 10; end %start ploting at tmin secs

rec_length = 120; %seconds to plot for each channel


par_macro.sr = 2000;
par_macro.detect_fmin = 1;
par_macro.detect_fmax = 120;
par_macro.auto = 1;

par_micro.sr = 30000;
par_micro.detect_fmin = 300;
par_micro.detect_fmax = 3000;
par_micro.auto = 0;

close all

for_text = 10;

load('NSx','NSx');

% poschs = find();

NSx = NSx(ismember(cell2mat({NSx.chan_ID}),channels));

bundles_to_plot = unique({NSx.bundle});
%     fig_num = ceil(10000*rand(1));
fig_num = 1032;
figure(fig_num)
set(fig_num, 'PaperUnits', 'inches', 'PaperType', 'A4', 'PaperPositionMode', 'auto','units','normalized','outerposition',[0 0 1 1], 'Visible', 'off') 

for ibun = 1:length(bundles_to_plot)
    pos_chans_to_plot = find(arrayfun(@(x) (strcmp(x.bundle,bundles_to_plot{ibun})),NSx));
    posch = pos_chans_to_plot(1);
        
    fprintf('%s ,', NSx(posch).bundle)

    % posch = find(arrayfun(@(x) (x.chan_ID==channels(1)),NSx));
    if NSx(posch).is_micro
        par = par_micro;
        mVmin = 50;
        w_pre=20;                       %number of pre-event data points stored
        w_post=44;                      %number of post-event data points stored
        min_ref_per=1.5;                                    %detector dead time (in ms)
        ref = floor(min_ref_per*par.sr/1000);                  %number of counts corresponding the dead time
        par.ref = ref;
        factor_thr=5;
    elseif NSx(posch).sr ==2000
        par = par_macro;
    end

    max_subplots = length(pos_chans_to_plot);
    Vmin = zeros(1,max_subplots);
    Vmax = zeros(1,max_subplots);
    channel_text_macro = cell(1,max_subplots);

    if NSx(posch).lts<par.sr * tmin
        disp('tmin is smaller than the recording length')
    else
        min_record = par.sr * tmin;
    end
    max_record = floor(min(NSx(posch).lts,min_record + par.sr * rec_length));
    tmax = max_record/par.sr;
    
    [b_orig,a_orig]=ellip(filt_order,0.1,40,[par.detect_fmin par.detect_fmax]*2/(par.sr));

    cont=1;    
    clf(fig_num)    
    for k= 1:length(pos_chans_to_plot)
        
        % LOAD NSX DATA
        channel1=NSx(pos_chans_to_plot(k)).chan_ID;
        posch = pos_chans_to_plot(k);
    
        if isfield(NSx,'dc') && ~isempty(NSx(posch).dc)
            dc = NSx(posch).dc;
        else
            dc=0;
        end
        f1 = fopen(sprintf('%s%s',NSx(posch).output_name,NSx(posch).ext),'r','l');
        fseek(f1,(min_record-1)*2,'bof');
        Samples = fread(f1,(max_record-min_record+1),'int16=>double')*NSx(posch).conversion + dc;
        fclose(f1);
    
        fprintf('%d,',channel1)
        b = b_orig;
        a = a_orig;
        if notchfilter
            [~, process_info] = pre_processing([],channel1);
            if ~isempty(process_info)
    %             [sos,g] = tf2sos(b,a);
                [sos,g] = tf2sos(b_orig,a_orig);
                g = g * process_info.G;
                sos = [process_info.SOS; sos];
                b = sos;
                a = g;
            end
        end
        % HIGH-PASS FILTER OF THE DATA
        xd=fast_filtfilt(b,a,Samples);
        
        clear Samples;
        
        if NSx(posch).is_micro
            %     % GET THRESHOLD AND NUMBER OF SPIKES BETWEEN 0 AND TMAX
            thr = factor_thr * median(abs(xd))/0.6745;
            thrmax = 10 * thr;     %thrmax for artifact removal is based on sorted settings.
            Vlim = mVmin;
        end
    %     xaux = find((abs(xd(w_pre+2:end-w_post-2)) > abs(thr)) & (abs(xd(w_pre+2:end-w_post-2)) < abs(thrmax))) +w_pre+1;
    %     nspk=nnz(diff(xaux)>ref)+1;
       
        % MAKES PLOT
        subplot(max_subplots,1,cont)
        box off; hold on        
        eje = linspace(tmin,tmax,length(xd));
        line(eje,xd)
        yl = ylim;
%         Vmin(cont) = yl(1);
%         Vmax(cont) = yl(2);

        Vmin(cont) = 1.05*prctile(xd,0.5);
        Vmax(cont) = 1.05*prctile(xd,99.5);
         
        ylabel(['Ch.' num2str(NSx(posch).chan_ID)])
        if NSx(posch).is_micro
            switch detect
                case 'pos'
                    xaux = find((xd(w_pre+2:end-w_post-2) > thr) & (abs(xd(w_pre+2:end-w_post-2)) < thrmax)) +w_pre+1;
                    line([tmin tmax],[thr thr],'color','r')
                    if ~par.auto
                        ylim([-Vlim 2*Vlim])
                    end
                case 'neg'
                    xaux = find((xd(w_pre+2:end-w_post-2) < -thr) & (abs(xd(w_pre+2:end-w_post-2)) < thrmax)) +w_pre+1;
                    line([tmin tmax],[-thr -thr],'color','r')
                    if ~par.auto
                        ylim([-2*Vlim Vlim])
                    end
                case 'both'
                    xaux = find((abs(xd(w_pre+2:end-w_post-2)) > thr) & (abs(xd(w_pre+2:end-w_post-2)) < abs(thrmax))) +w_pre+1;
                    line([tmin tmax],[thr thr],'color','r')
                    line([tmin tmax],[-thr -thr],'color','r')
                    if ~par.auto
                        ylim([-2*Vlim 2*Vlim])
                    end
            end
            nspk=nnz(diff(xaux)>ref)+1;
            if (nspk > ceil(rec_length/20) && nspk < rec_length * 60)
                ylabel(['Ch.' num2str(NSx(posch).chan_ID)],'fontsize',10,'fontweight','bold')
            end
            if notchfilter==0  || nnz(arrayfun(@(x) (x.chID==channel1),process_info))==0
                text((tmax-tmin)/2+tmin+5,for_text+max(ylim),sprintf('%s. %d spikes',NSx(posch).label,nspk),'fontsize',10)
            elseif nnz(arrayfun(@(x) (x.chID==channel1),process_info))>0
                posch_notch = find(arrayfun(@(x) (x.chID==channel1),process_info));
                if isempty(posch_notch)
                    cant_notches1l = 0; cant_notches1h = 0;
                else
                    cant_notches1l = sum(process_info(posch_notch).freqs_notch<par.detect_fmin);
                    cant_notches1h = sum(process_info(posch_notch).freqs_notch>par.detect_fmin & process_info(posch_notch).freqs_notch<par.detect_fmax);
                end
                text((tmax-tmin)/3+tmin,for_text+max(ylim),sprintf('%s. %d spikes. %d and %d notches applied below and above %d Hz',NSx(posch).label,nspk,cant_notches1l,cant_notches1h,par.detect_fmin),'fontsize',10)
            end
        else
            posch_notch = find(arrayfun(@(x) (x.chID==channel1),process_info));
            if isempty(posch_notch)
                cant_notches1l = 0; cant_notches1h = 0; 
            else
                cant_notches1l = sum(process_info(posch_notch).freqs_notch<200);
                cant_notches1h = sum(process_info(posch_notch).freqs_notch>200 & process_info(posch_notch).freqs_notch<1000);
            end
            channel_text_macro{cont} = sprintf('%s. %d and %d notches applied below and above %d Hz',NSx(posch).label,cant_notches1l,cant_notches1h,200);
        end
        clear xd;         
        
        set(gca,'Xlim',[tmin tmax],'Xtick',linspace(tmin,tmax,7))
                
        if cont==max_subplots
            fprintf('\n')
            xlabel('Time (sec)')
            bundle = NSx(posch).bundle;

            if ~NSx(posch).is_micro && par.auto
                for jj=1:max_subplots
                    subplot(max_subplots,1,jj)
                    ylim([min(Vmin) max(Vmax)])
%                     text((tmax-tmin)/3+tmin,for_text+max(ylim),channel_text_macro{jj},'fontsize',10)
                    text((tmax-tmin)/3+tmin,0.9*15*max_subplots*max(ylim)/diff(ylim)/4+max(ylim),channel_text_macro{jj},'fontsize',10)
                end
            end

            currentFolder = pwd;
            [~, folderName] = fileparts(currentFolder);                    
            title_out = sprintf('%s.   bundle  %s.   fmin %d Hz. fmax %d Hz',folderName,bundle,par.detect_fmin,par.detect_fmax);
            
            sgtitle(title_out,'fontsize',12,'interpreter','none','fontWeight','bold','HorizontalAlignment','center')
            
            if NSx(posch).is_micro
                outfile=sprintf('%s_%s_filtorder%d_withnotches%d_det%s','fig2print_bundle',bundle,filt_order,notchfilter,detect);
            else
                outfile=sprintf('%s_%s_filtorder%d_withnotches%d','fig2print_bundle',bundle,filt_order,notchfilter);
            end
            print(fig_num,fullfile(pwd,outfile),'-dpng')
            
            clf
        end 

        cont=cont+1;
    end
end
close(fig_num)

    