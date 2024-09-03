function save_grating_stimuli_powerspectrum_folder(curfolder)
% for all grating_stimuli get the part of the trace
% to pd & nd directions
param_speed = 5;
param_spf = 0.01;
Fs = 10000; %sampling rate
stim_dur = 1*Fs; 
%exclude motion onset time
nexcl = 0.1*Fs;
%approximate duration of the traces
L_appr=stim_dur-nexcl; 
% frequency bins
f =  fft_bins(Fs,L_appr);  
% we don't care about frequencies higher than 1K
fmax=10^3;
fmaxind=find(f<=fmax,1,'last');
f=f(1:fmaxind);

% curfolder='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpND\VS\VS1_4\210128';
    %% get strain, cell type and date from the folder path
    [strain,cellgroup,celltype,datestr,cellid] = cell_info_from_path(curfolder);
    
    %% get all stimulus files
    filelist = dir(curfolder);%get list of files
    filelist = filelist(~[filelist.isdir]);  %remove folders from list

    stim_name='_grating_stimuli_';
    
    %go thru the list and collect files with'_grating_stimuli_' 
    stim_files={};
    i=1;
    for li =1:length(filelist)
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind) && ...
                contains(filelist(li).name,'.log','IgnoreCase', true) && ...
                ~(contains(filelist(li).name, '_Vel') || ... %exclude velocity tuning curves
                  contains(filelist(li).name, '_expansion_') || ... %exclude grating expansion
                  contains(filelist(li).name, '_contrast_levels_') ||... %exclude grating contrast levels
                  contains(filelist(li).name, '_split_screen_')) %exclude grating split screen
            stim_files{i}=filelist(li).name;            
            i=i+1;
        end
    end
    
    resfolder=fullfile(curfolder,'res');
    if ~exist(resfolder,'dir')
        mkdir(resfolder);
    end
    
    total_nrep=0;
    %% iterate thru stim files, get raw PD/ND traces, make spectrum and average
    nfiles=length(stim_files);
    ncount=0;
    logiles={};
    prfiles={};
    for i=1:nfiles        
        fulllog_file=fullfile(curfolder,stim_files{i});
        stimparam_i= load_grating_stimparam(fulllog_file);
        if stimparam_i.sp_freq~=param_spf || stimparam_i.speed~=param_speed || stimparam_i.baseline_dur<1
            continue;
        end
                
        rawdata_i=trace_pd_nd_grating_one_recording(fulllog_file); 
        if isempty(rawdata_i)            
            continue;
        end
        ncount=ncount+1;
        %structure to pool all data        
        logiles{ncount}=stim_files{i};        
        prfiles{ncount}=rawdata_i.prname;
        data(ncount).trace=rawdata_i.trace;        
        nreps=size(data(ncount).trace.pd_traces,1);
        total_nrep=total_nrep+nreps;
    end

    if isempty(logiles)
        return;
    end
    
    %% fft power spectrum for all traces
    powerspectrum_pd=zeros(total_nrep,fmaxind);    
    powerspectrum_nd=zeros(total_nrep,fmaxind);    
    powerspectrum_blpd_a=zeros(total_nrep,fmaxind);
    powerspectrum_blpd_aa=zeros(total_nrep,fmaxind);
    powerspectrum_blnd_a=zeros(total_nrep,fmaxind);
    powerspectrum_blnd_aa=zeros(total_nrep,fmaxind);
    
    %contribution of parts of spectrum: low, mid, total
    pow_pd = zeros(total_nrep,3);
    pow_nd = zeros(total_nrep,3);
    pow_blpd_a = zeros(total_nrep,3);
    pow_blpd_aa = zeros(total_nrep,3);
    pow_blnd_a = zeros(total_nrep,3);
    pow_blnd_aa = zeros(total_nrep,3);
    
    %% for each rep and recording make a powerspectrum and average
    figure, t = tiledlayout(3,1,'TileSpacing','compact');
    ax1=nexttile; hold(ax1); 
    ax2=nexttile; hold(ax2); 
    ax3=nexttile; hold(ax3); 
    nall=1;    
    for i=1:ncount
        nrepi = size(data(i).trace.pd_traces,1);
        for re=1:nrepi
            nlen=min(stim_dur,size(data(i).trace.pd_traces(re,:),2));
            pd_reptrace = data(i).trace.pd_traces(re,1:nlen);
            plot(ax1,pd_reptrace);
            pd_reptrace = pd_reptrace(nexcl+1:end);
            pd_reptrace=pd_reptrace-mean(pd_reptrace);

            nlen=min(stim_dur,size(data(i).trace.nd_traces(re,:),2));
            nd_reptrace = data(i).trace.nd_traces(re,1:nlen);
            plot(ax2,nd_reptrace);
            nd_reptrace = nd_reptrace(nexcl+1:end);
            nd_reptrace=nd_reptrace-mean(nd_reptrace);

            bl_pd_reptrace = data(i).trace.pd_baseline(re,:);
            plot(ax3,bl_pd_reptrace);
            bl_pd_reptrace=bl_pd_reptrace-mean(bl_pd_reptrace);
            bl_nd_reptrace = data(i).trace.nd_baseline(re,:);
            plot(ax3,bl_nd_reptrace);
            bl_nd_reptrace=bl_nd_reptrace-mean(bl_nd_reptrace);
            
            [powerspectrum_pd(nall,:), pow_pd(nall,1), pow_pd(nall,2), pow_pd(nall,3),~] = get_trace_pwspectrum_maxlen(pd_reptrace, Fs, L_appr, fmax);
            [powerspectrum_nd(nall,:), pow_nd(nall,1), pow_nd(nall,2), pow_nd(nall,3),~] = get_trace_pwspectrum_maxlen(nd_reptrace, Fs, L_appr, fmax);
            blhalf=round(length(bl_pd_reptrace)/2);
            [powerspectrum_blpd_a(nall,:), pow_blpd_a(nall,1), pow_blpd_a(nall,2), pow_blpd_a(nall,3),~] = get_trace_pwspectrum_maxlen(bl_pd_reptrace(1:blhalf), Fs, L_appr, fmax);
            [powerspectrum_blpd_aa(nall,:), pow_blpd_aa(nall,1), pow_blpd_aa(nall,2), pow_blpd_aa(nall,3),~] = get_trace_pwspectrum_maxlen(bl_pd_reptrace(blhalf+1:end), Fs, L_appr, fmax);
            blhalf=round(length(bl_nd_reptrace)/2);
            [powerspectrum_blnd_a(nall,:), pow_blnd_a(nall,1), pow_blnd_a(nall,2), pow_blnd_a(nall,3),~] = get_trace_pwspectrum_maxlen(bl_nd_reptrace(1:blhalf), Fs, L_appr, fmax);
            [powerspectrum_blnd_aa(nall,:), pow_blnd_aa(nall,1), pow_blnd_aa(nall,2), pow_blnd_aa(nall,3),~] = get_trace_pwspectrum_maxlen(bl_nd_reptrace(blhalf+1:end), Fs, L_appr, fmax);
            nall=nall+1;            
        end
    end

    ylabel(ax1,'PD');
    ylabel(ax2,'ND');
    ylabel(ax3,'bl');
    title(t,strjoin({strain,cellgroup,celltype,datestr,cellid},' '),'Interpreter','none');

    %save averages and std    
    exp_data.pw.powerspectrum_pd_av=mean(powerspectrum_pd,1);
    exp_data.pw.powerspectrum_nd_av=mean(powerspectrum_nd,1);
    exp_data.pw.powerspectrum_blpd_a_av=mean(powerspectrum_blpd_a,1);  
    exp_data.pw.powerspectrum_blpd_aa_av=mean(powerspectrum_blpd_aa,1);  
    exp_data.pw.powerspectrum_blnd_a_av=mean(powerspectrum_blnd_a,1);  
    exp_data.pw.powerspectrum_blnd_aa_av=mean(powerspectrum_blnd_aa,1);  
    exp_data.pw.powerspectrum_pd_std=std(powerspectrum_pd,[],1);
    exp_data.pw.powerspectrum_nd_std=std(powerspectrum_nd,[],1);
    exp_data.pw.powerspectrum_blpd_a_std=std(powerspectrum_blpd_a,[],1);  
    exp_data.pw.powerspectrum_blpd_aa_std=std(powerspectrum_blpd_aa,[],1);  
    exp_data.pw.powerspectrum_blnd_a_std=std(powerspectrum_blnd_a,[],1);  
    exp_data.pw.powerspectrum_blnd_aa_std=std(powerspectrum_blnd_aa,[],1);  

    exp_data.pw.pow_pd_av = mean(pow_pd,1);    
    exp_data.pw.pow_pd_std = std(pow_pd,[],1);
    exp_data.pw.pow_nd_av = mean(pow_nd,1);    
    exp_data.pw.pow_nd_std = std(pow_nd,[],1);
    exp_data.pw.pow_blpd_a_av = mean(pow_blpd_a,1);    
    exp_data.pw.pow_blpd_a_std = std(pow_blpd_a,[],1);
    exp_data.pw.pow_blpd_aa_av = mean(pow_blpd_aa,1);    
    exp_data.pw.pow_blpd_aa_std = std(pow_blpd_aa,[],1);
    exp_data.pw.pow_blnd_a_av = mean(pow_blnd_a,1);    
    exp_data.pw.pow_blnd_a_std = std(pow_blnd_a,[],1);
    exp_data.pw.pow_blnd_aa_av = mean(pow_blnd_aa,1);    
    exp_data.pw.pow_blnd_aa_std = std(pow_blnd_aa,[],1);

    %normalized mean traces and range contributions    
    exp_data.pw_normed.powerspectrum_pd_av_normed=exp_data.pw.powerspectrum_pd_av/sum(exp_data.pw.powerspectrum_pd_av);
    exp_data.pw_normed.powerspectrum_nd_av_normed=exp_data.pw.powerspectrum_nd_av/sum(exp_data.pw.powerspectrum_nd_av);
    exp_data.pw_normed.powerspectrum_blpd_a_av_normed=exp_data.pw.powerspectrum_blpd_a_av/sum(exp_data.pw.powerspectrum_blpd_a_av);  
    exp_data.pw_normed.powerspectrum_blpd_aa_av_normed=exp_data.pw.powerspectrum_blpd_aa_av/sum(exp_data.pw.powerspectrum_blpd_aa_av);  
    exp_data.pw_normed.powerspectrum_blnd_a_av_normed=exp_data.pw.powerspectrum_blnd_a_av/sum(exp_data.pw.powerspectrum_blnd_a_av);  
    exp_data.pw_normed.powerspectrum_blnd_aa_av_normed=exp_data.pw.powerspectrum_blnd_aa_av/sum(exp_data.pw.powerspectrum_blnd_aa_av);

    exp_data.pw_normed.pow_pd_av = exp_data.pw.pow_pd_av/exp_data.pw.pow_pd_av(3);
    exp_data.pw_normed.pow_nd_av = exp_data.pw.pow_nd_av/exp_data.pw.pow_nd_av(3);
    exp_data.pw_normed.pow_blpd_a_av = exp_data.pw.pow_blpd_a_av/exp_data.pw.pow_blpd_a_av(3);
    exp_data.pw_normed.pow_blpd_aa_av = exp_data.pw.pow_blpd_aa_av/exp_data.pw.pow_blpd_aa_av(3);
    exp_data.pw_normed.pow_blnd_a_av = exp_data.pw.pow_blnd_a_av/exp_data.pw.pow_blnd_a_av(3);
    exp_data.pw_normed.pow_blnd_aa_av = exp_data.pw.pow_blnd_aa_av/exp_data.pw.pow_blnd_aa_av(3);

    exp_data.pw.nrep= total_nrep;
    exp_data.pw.freq=f;
    
    

                
%     figure;
%     loglog(f,powerspectrum_av);
    %power log to stretch the low range
%     plot(log10(f),log10(powerspectrum_av));
%     flog=log10(f);    
%     ftick=linspace(flog(2),flog(end),5);
%     XTickLabels = cellstr(num2str(round(ftick(:)), '10^{%d}'));
%     xticks(ftick);
%     xticklabels(XTickLabels);
%     ylog=[min(log10(psp)),max(log10(psp))];   
%     ytick=linspace(ylog(1),ylog(2),5);
%     YTickLabels = cellstr(num2str(round(ytick(:)), '10^{%d}'));
%     yticks(ytick);
%     yticklabels(YTickLabels);
    
   

    
    %create data structure    
    exp_data.folder=curfolder;
    exp_data.logfiles=logiles;
    exp_data.prfiles=prfiles;    
    %cell info
    exp_data.cellinfo.strain = strain;
    exp_data.cellinfo.cellgroup = cellgroup;
    exp_data.cellinfo.celltype = celltype;
    exp_data.cellinfo.datestr = datestr;
    exp_data.cellinfo.cellid = cellid;       

    %save    
    if isempty(cellid)
        filename_i=['grating_pd_nd_pwspectrum_',strain,'_',celltype,'_',datestr];        
    else
        filename_i=['grating_pd_nd_pwspectrum_',strain,'_',celltype,'_',datestr,'_',cellid];        
    end    
    fullfile_i=fullfile(resfolder,[filename_i,'.mat']);
    save(fullfile_i,'exp_data','-v7.3');  
end


% 
% 
% figure, 
%     for ri=1:nrep
%         plot(trace_all(ri).reptrace);
%         hold on;
%     end
%     
%     figure, plot(reptrace);
% 
%     Fs = 10000;          %sampling rate
%     T = 1/Fs;           %sampling interval
%     L = numel(reptrace);            %Number of time points
%     t = 0:T:(L-1)*T;    %time vector
%     
% 
%     figure;
%     plot(t*1000,reptrace)
%     xlabel('time (in ms)')
%     ylabel('mV')
%     title('Original Signal')
% 
%     N = Fs/2*1000;%4096;%2048;%1024;           %FFT points
%     % FFT
%     X = fft(reptrace,N);
%     %make one sided signal
%     SSB = X(1:N/2);
%     SSB(2:end) = 2*SSB(2:end);
%     %convert bins to frequences
%     f = (0:N/2-1)*(Fs/N);
%     % Amplitude
% %     figure;
% %     plot(f,abs(SSB/L))
% %     psp = abs(SSB).^2/(L^2);
% 
%     psp = (abs(SSB)/L*Fs).^2; 
% %     %power is the square of the amplitude of frequencies
% %     psp = (abs(SSB)/L).^2;
% %     %in borst paper the units are V^2s
% %     psp = ((abs(SSB)/L).^2)*Fs;
% 
%     figure;
%     %power log to stretch the low range
%     plot(log10(f),log10(psp));
%     flog=log10(f);    
%     ftick=linspace(flog(2),flog(end),5);
%     XTickLabels = cellstr(num2str(round(ftick(:)), '10^{%d}'));
%     xticks(ftick);
%     xticklabels(XTickLabels);
%     ylog=[min(log10(psp)),max(log10(psp))];   
%     ytick=linspace(ylog(1),ylog(2),5);
%     YTickLabels = cellstr(num2str(round(ytick(:)), '10^{%d}'));
%     yticks(ytick);
%     yticklabels(YTickLabels);
%     title('FlpND')
%     
% %     plot(f(1:idcut),abs(SSB(1:idcut)/L))
% % xlabel('f (in Hz)')
% % ylabel('|X(f)|')
% 
% 
% 
% 
%     figure, 
%     [s,f,t,ps] = spectrogram(reptrace,1024,512,[0:1:250],10000,'yaxis');
%     s0_10 = sum(s(1:10,:),1);
%     figure, plot(t, s0_10)
%     title('Power per 0-10 freq');
% 
%     s10_20 = sum(s(11:20,:),1);
%     figure, plot(t, s10_20)
%     title('Power per 10-20 freq');
% 
%     psn=ps/sum(ps(:));
%     sum(psn(1:10,:),'all')

function [psp, pow_lowfft, pow_midfft, pow_totalfft] = get_trace_amplitude_spectrum(reptrace, sampling_rate)        
    L = numel(reptrace); %Number of time points    
    N = 2^nextpow2(L); %sampling_rate/2*1000;%4096;   %FFT points
    f =  fft_bins(sampling_rate, L); 
    %low freq range [0,10];
    f10=find(f<=10,1,'last');
    f50=find(f<=50,1,'last');
    
    % FFT
    X = fft(reptrace,N);
    %make one sided signal
    SSB = X(1:N/2);
    SSB(2:end) = 2*SSB(2:end);            
    % Amplitude
    psp = abs(SSB/L); 

    %total power within low freq <10Hz    
    pow_lowfft=sum(psp(1:f10));
    %total power within mid freq [10Hz, 50Hz]
    pow_midfft=sum(psp(f10+1:f50));
    pow_totalfft = sum(psp);        
end

function [psp, pow_lowfft, pow_midfft, pow_totalfft, f] = get_trace_pwspectrum_maxlen(reptrace, sampling_rate, Lmax, maxfreq)        
    L = numel(reptrace); %length of the signal 
    %N = sampling_rate/2*1000;%2^nextpow2(L);%4096;%2048;%1024;      %FFT points
    N = 2^nextpow2(Lmax);
    f =  fft_bins(sampling_rate, Lmax); 
    %low freq range [0,10];
    f10=find(f<=10,1,'last');
    f50=find(f<=50,1,'last');
    if ~exist('maxfreq','var') || isempty(maxfreq)
        maxfreq=Inf;
        fmaxind=length(f);
    else
        fmaxind = find(f<=maxfreq,1,'last');
    end
               
    % FFT 
    X = fft(reptrace,N);
    %take one side only
    X = X(1:N/2+1);

    % power normalized by freq resolution
    psp = (1/(sampling_rate*L)) * abs(X).^2;
    % adjust because one-sided
    psp(2:end-1) = 2*psp(2:end-1);

    %total power within low freq <10Hz    
    pow_lowfft=sum(psp(1:f10));
    %total power within mid freq [10Hz, 50Hz]
    pow_midfft=sum(psp(f10+1:f50));
    pow_totalfft = sum(psp); 

    %crop frequences to the max specified
    f=f(1:fmaxind);
    psp=psp(1:fmaxind);
end

% function [psp, pow_lowfft, pow_midfft, pow_totalfft, f] = get_trace_pwspectrum(reptrace, sampling_rate, maxfreq)        
%     L = numel(reptrace); %length of the signal 
%     %N = sampling_rate/2*1000;%2^nextpow2(L);%4096;%2048;%1024;      %FFT points
%     N = 2^nextpow2(L);
%     f =  fft_bins(sampling_rate, L); 
%     %low freq range [0,10];
%     f10=find(f<=10,1,'last');
%     f50=find(f<=50,1,'last');
%     if ~exist('maxfreq','var') || isempty(maxfreq)
%         maxfreq=Inf;
%         fmaxind=length(f);
%     else
%         fmaxind = find(f<=maxfreq,1,'last');
%     end
%                
%     % FFT 
%     X = fft(reptrace,N);
%     %take one side only
%     X = X(1:N/2+1);
% 
%     % power normalized by freq resolution
%     psp = (1/(sampling_rate*L)) * abs(X).^2;
%     % adjust because one-sided
%     psp(2:end-1) = 2*psp(2:end-1);
% 
%     %total power within low freq <10Hz    
%     pow_lowfft=sum(psp(1:f10));
%     %total power within mid freq [10Hz, 50Hz]
%     pow_midfft=sum(psp(f10+1:f50));
%     pow_totalfft = sum(psp); 
% 
%     %crop frequences to the max specified
%     f=f(1:fmaxind);
%     psp=psp(1:fmaxind);
% end

function f = fft_bins(Fs, L)  
    N=2^nextpow2(L);
%     N = Fs/2*1000;%2^nextpow2(L);
    %convert FFT bins to frequences
    f = 0:(Fs/N):(Fs/2);
end
