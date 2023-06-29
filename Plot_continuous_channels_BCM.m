function Plot_continuous_channels_BCM(channels,notchfilter,filt_order,detect)
% Plots the continuous data and does a first detection of spikes for each
% channel. 
detname = 1; 
if ~exist('filt_order','var')|| isempty(filt_order),     filt_order = 4; end
if ~exist('notchfilter','var')|| isempty(notchfilter),     notchfilter = 0; end
if ~exist('detect','var')|| isempty(detect),     detect = 'both'; detname = 0; end  %'pos','neg','both'

auto = 0;
tmin = 10;%start ploting at tmin secs
rec_length = 120; %seconds to plot for each channel
par.sr = 30000;
par.detect_fmin = 300;
par.detect_fmax = 3000;

close all
Vmin = 50;
for_text = 10;
w_pre=20;                       %number of pre-event data points stored
w_post=44;                      %number of post-event data points stored
min_ref_per=1.5;                                    %detector dead time (in ms)
ref = floor(min_ref_per*par.sr/1000);                  %number of counts corresponding the dead time
par.ref = ref;
factor_thr=5;

load('NSx','NSx');
posch = find(arrayfun(@(x) (x.chan_ID==channels(1)),NSx));

if NSx(posch).lts<par.sr * tmin
    disp('tmin is smaller than the recording length')
else
    min_record = par.sr * tmin;
end
max_record = floor(min(NSx(posch).lts,min_record + par.sr * rec_length));
tmax = max_record/par.sr;

max_subplots = 8;

fig_num = ceil(10000*rand(1));
figure(fig_num)
set(fig_num, 'PaperUnits', 'inches', 'PaperType', 'A4', 'PaperPositionMode', 'auto','units','normalized','outerposition',[0 0 1 1], 'Visible', 'off') 

cont=1;

clf(fig_num)
[b_orig,a_orig]=ellip(filt_order,0.1,40,[par.detect_fmin par.detect_fmax]*2/(par.sr));

ch = [];

for k= 1:length(channels)
    % LOAD NSX DATA
    channel1=channels(k);
    ch = [ch channel1];
    posch = find(arrayfun(@(x) (x.chan_ID==channel1),NSx));
    if isfield(NSx,'dc') && ~isempty(NSx(posch).dc)
        dc = NSx(posch).dc;
    else
        dc=0;
    end
    f1 = fopen(sprintf('%s%s',NSx(posch).output_name,NSx(posch).ext),'r','l');
    fseek(f1,(min_record-1)*2,'bof');
    Samples = fread(f1,(max_record-min_record+1),'int16=>double')*NSx(posch).conversion+dc;
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
    
    % GET THRESHOLD AND NUMBER OF SPIKES BETWEEN 0 AND TMAX
    thr = factor_thr * median(abs(xd))/0.6745;
    thrmax = 10 * thr;     %thrmax for artifact removal is based on sorted settings.
    
%     xaux = find((abs(xd(w_pre+2:end-w_post-2)) > abs(thr)) & (abs(xd(w_pre+2:end-w_post-2)) < abs(thrmax))) +w_pre+1;
%     nspk=nnz(diff(xaux)>ref)+1;
   
    % MAKES PLOT
    subplot(max_subplots,1,cont)
    box off; hold on
    Vlim = Vmin;
    eje = linspace(tmin,tmax,length(xd));
    line(eje,xd)
    switch detect
        case 'pos'
            xaux = find((xd(w_pre+2:end-w_post-2) > thr) & (abs(xd(w_pre+2:end-w_post-2)) < thrmax)) +w_pre+1;
            line([tmin tmax],[thr thr],'color','r')
            if ~auto
                ylim([-Vlim 2*Vlim])
            end
        case 'neg'
            xaux = find((xd(w_pre+2:end-w_post-2) < -thr) & (abs(xd(w_pre+2:end-w_post-2)) < thrmax)) +w_pre+1;
            line([tmin tmax],[-thr -thr],'color','r')
            if ~auto
                ylim([-2*Vlim Vlim])
            end
        case 'both'
            xaux = find((abs(xd(w_pre+2:end-w_post-2)) > thr) & (abs(xd(w_pre+2:end-w_post-2)) < abs(thrmax))) +w_pre+1;
            line([tmin tmax],[thr thr],'color','r')
            line([tmin tmax],[-thr -thr],'color','r')
            if ~auto
                ylim([-2*Vlim 2*Vlim])
            end
    end
    clear xd;    
    nspk=nnz(diff(xaux)>ref)+1;

    ylabel(['Ch.' num2str(NSx(posch).chan_ID)])        

    cont=cont+1;
    if (nspk > ceil(rec_length/20) && nspk < rec_length * 60)
        ylabel(['Ch.' num2str(NSx(posch).chan_ID)],'fontsize',10,'fontweight','bold')
    end
    if notchfilter==0  || nnz(arrayfun(@(x) (x.chID==channel1),process_info))==0
        text((tmax-tmin)/2+tmin+5,for_text+max(ylim),[num2str(nspk) '  spikes.'],'fontsize',10)
    elseif nnz(arrayfun(@(x) (x.chID==channel1),process_info))>0
        posch_notch = find(arrayfun(@(x) (x.chID==channel1),process_info));
        cant_notches1l = sum(process_info(posch_notch).freqs_notch<par.detect_fmin);
        cant_notches1h = sum(process_info(posch_notch).freqs_notch>par.detect_fmin & process_info(posch_notch).freqs_notch<par.detect_fmax);
        text((tmax-tmin)/3+tmin,for_text+max(ylim),[num2str(nspk) '  spikes. ' num2str(cant_notches1l) ' and ' num2str(cant_notches1h) ' notches applied below and above ' num2str(par.detect_fmin) ' Hz'],'fontsize',10)
    end
    set(gca,'Xlim',[tmin tmax],'Xtick',linspace(tmin,tmax,7))
    if(~mod(cont-1,max_subplots))
        fprintf('\n')
        cont=1;
        xlabel('Time (sec)')
        macro = NSx(posch).macro;
        title_out = sprintf('%s.   bundle  %s.   fmin %d Hz. fmax %d Hz',pwd,macro,par.detect_fmin,par.detect_fmax);
        
        sgtitle(title_out,'fontsize',12,'interpreter','none','fontWeight','bold','HorizontalAlignment','left')
        
        if detname
            outfile=sprintf('%s_%s_filtorder%d_withnotches%d_det%s','fig2print_macro',macro,filt_order,notchfilter,detect);
        else
            outfile=sprintf('%s_%s_filtorder%d_withnotches%d','fig2print_macro',macro,filt_order,notchfilter);
        end
        print(fig_num,fullfile(pwd,outfile),'-dpng')
        
        
        ch = [];
        clf
    end 
end
close(fig_num)

    