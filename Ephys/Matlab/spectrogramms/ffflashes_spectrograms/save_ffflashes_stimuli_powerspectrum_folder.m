function save_ffflashes_stimuli_powerspectrum_folder(curfolder)
% for all full field flashes folder get the part of the traces
% correspondign to on and off flash
Fs = 10000; %sampling rate
flash_dur_s=2;
tbefore_s=1;
tafter_s=1;

stim_dur = flash_dur_s*Fs; 
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

    stim_name='_full_field_flashes_';
    
    %go thru the list and collect files with'_full_field_flashes_' 
    stim_files={};
    i=1;
    for li =1:length(filelist)
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind)
            stim_files{i}=filelist(li).name;            
            i=i+1;
        end
    end
    
    resfolder=fullfile(curfolder,'res');
    if ~exist(resfolder,'dir')
        mkdir(resfolder);
    end
    
    total_nrep=0;
    total_nrep_off=0;
    %% iterate thru stim files, get raw ON/OFF traces, make spectrum and average
    nfiles=length(stim_files);
    ncount=0;
    logiles={};
    prfiles={};
    for i=1:nfiles        
        fulllog_file=fullfile(curfolder,stim_files{i});
                
        rawdata_i=trace_on_of_flashes_one_recording(fulllog_file, flash_dur_s, tbefore_s, tafter_s); 
        if isempty(rawdata_i)            
            continue;
        end
        ncount=ncount+1;
        %structure to pool all data        
        logiles{ncount}=stim_files{i};        
        prfiles{ncount}=rawdata_i.prname;
        data(ncount).trace=rawdata_i.trace;                
        total_nrep=total_nrep+size(data(ncount).trace.on_traces,1);
        total_nrep_off=total_nrep_off+size(data(ncount).trace.off_traces,1);
    end

    if isempty(logiles)
        return;
    end    
    
    %%to check the data plot raw extracted traces during and off flash
    figure;
    ax1 = subplot(3,4,1:4);
    ax2 = subplot(3,4,5:8);
    ax3 = subplot(3,4,9:12);

    for i=1:ncount
        nrepi = size(data(i).trace.on_traces,1);
        hold(ax1);       
        for re=1:nrepi            
            on_reptrace = data(i).trace.on_traces(re,:);
            plot(ax1,on_reptrace); 
        end
        hold(ax1);       
        nrepi = size(data(i).trace.off_traces,1);
        hold(ax2);       
        for re=1:nrepi            
            off_reptrace = data(i).trace.off_traces(re,:);
            plot(ax2,off_reptrace);  
        end 
        hold(ax2);  
        hold(ax3);       
        for re=1:nrepi            
            beforeflash_reptrace = data(i).trace.beforeflash_traces(re,:);
            plot(ax3,beforeflash_reptrace); 
        end
        hold(ax3);       
    end
    ylabel(ax1,'ON');
    ylabel(ax2,'OFF');
    ylabel(ax3,'before ON flash');
    sgtitle(curfolder,'Interpreter','none');
    linkaxes([ax1,ax2],'x');


    %% for each rep and recording make a powerspectrum and average
    %% fft power spectrum for all traces    
    powerspectrum_on=zeros(total_nrep,fmaxind);    
    powerspectrum_off=zeros(total_nrep_off,fmaxind);   
    powerspectrum_before_on=zeros(total_nrep_off,fmaxind);   
    
    %contribution of parts of spectrum: low, mid, total
    pow_on = zeros(total_nrep,3);
    pow_off = zeros(total_nrep_off,3);
    pow_before_on = zeros(total_nrep_off,3);

    nall_on=1;
    nall_off=1;
    for i=1:ncount
        nrepi = size(data(i).trace.on_traces,1);
        for re=1:nrepi            
            %power spectrum of the ON-trace
            on_reptrace = data(i).trace.on_traces(re,:);
            on_reptrace = on_reptrace(nexcl+1:end);
            [powerspectrum_on(nall_on,:), pow_on(nall_on,1), pow_on(nall_on,2), pow_on(nall_on,3),~] = get_trace_pwspectrum_maxlen(on_reptrace, Fs, L_appr, fmax);            

            %power spectrum of the trace before on flash
            before_on_reptrace = data(i).trace.beforeflash_traces(re,:);            
            [powerspectrum_before_on(nall_on,:), pow_before_on(nall_on,1), pow_before_on(nall_on,2), pow_before_on(nall_on,3),~] = get_trace_pwspectrum_maxlen(before_on_reptrace, Fs, L_appr, fmax);
            nall_on=nall_on+1;
        end

        nrepi = size(data(i).trace.off_traces,1);
        for re=1:nrepi            
            off_reptrace = data(i).trace.off_traces(re,:);
            off_reptrace = off_reptrace(nexcl+1:end);
            [powerspectrum_off(nall_off,:), pow_off(nall_off,1), pow_off(nall_off,2), pow_off(nall_off,3),~] = get_trace_pwspectrum_maxlen(off_reptrace, Fs, L_appr, fmax);
            nall_off=nall_off+1;
        end 
    end

    %save averages and std    
    exp_data.pw.powerspectrum_on_av=mean(powerspectrum_on,1);
    exp_data.pw.powerspectrum_off_av=mean(powerspectrum_off,1);
    exp_data.pw.powerspectrum_before_on_av=mean(powerspectrum_before_on,1);
    exp_data.pw.powerspectrum_on_std=std(powerspectrum_on,[],1);
    exp_data.pw.powerspectrum_off_std=std(powerspectrum_off,[],1);
    exp_data.pw.powerspectrum_before_on_std=std(powerspectrum_before_on,[],1);    
    
    exp_data.pw.pow_on_av = mean(pow_on,1);    
    exp_data.pw.pow_on_std = std(pow_on,[],1);
    exp_data.pw.pow_off_av = mean(pow_off,1);    
    exp_data.pw.pow_off_std = std(pow_off,[],1);
    exp_data.pw.pow_before_on_av = mean(pow_before_on,1);    
    exp_data.pw.pow_before_on_std = std(pow_before_on,[],1);

    %normalized mean traces and range contributions    
    exp_data.pw_normed.powerspectrum_on_av_normed=exp_data.pw.powerspectrum_on_av/sum(exp_data.pw.powerspectrum_on_av);
    exp_data.pw_normed.powerspectrum_off_av=exp_data.pw.powerspectrum_off_av/sum(exp_data.pw.powerspectrum_off_av);
    exp_data.pw_normed.powerspectrum_before_on_av=exp_data.pw.powerspectrum_before_on_av/sum(exp_data.pw.powerspectrum_before_on_av);    
    
    exp_data.pw_normed.pow_on_av = exp_data.pw.pow_on_av/exp_data.pw.pow_on_av(3);  
    exp_data.pw_normed.pow_off_av = exp_data.pw.pow_off_av/exp_data.pw.pow_off_av(3);  
    exp_data.pw_normed.pow_before_on_av = exp_data.pw.pow_before_on_av/exp_data.pw.pow_before_on_av(3);      
    
    
    exp_data.pw.nrep_on= total_nrep;
    exp_data.pw.nrep_off= total_nrep_off;
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
        filename_i=['flashes_on_off_pwspectrum_',strain,'_',celltype,'_',datestr];        
    else
        filename_i=['flashes_on_off_pwspectrum_',strain,'_',celltype,'_',datestr,'_',cellid];        
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

%     f1000=find(f<=1000,1,'last');
%     figure, loglog(f(1:f1000),psp(1:f1000));
%     title('Pw');
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

function [psp,plow,pmid,ptot,f] =  pwspectrum_prev_style(trace)
    Fs = 10000;          %sampling rate
    T = 1/Fs;           %sampling interval
    N = Fs/2*1000;%4096;%2048;%1024;      %FFT points
    %convert FFT bins to frequences
    f = (0:N/2-1)*(Fs/N);
    f10=find(f<=10,1,'last');
    f50=find(f<=50,1,'last');
        
    %% make a powerspectrum
    L = numel(trace); %Number of time points
    
    % FFT
    X = fft(trace,N);
    %make one sided signal
    SSB = X(1:N/2);
    SSB(2:end) = 2*SSB(2:end);            
    % Amplitude
    psp = abs(SSB/L);
    plow=sum(psp(1:f10));
    pmid=sum(psp(10:f50));
    ptot=sum(psp);
    f1000=find(f<=1000,1,'last');

    figure, loglog(f(1:f1000),psp(1:f1000));
    title('Amplitude');
end


function f = fft_bins(Fs, L)  
    N=2^nextpow2(L);
%     N = Fs/2*1000;%2^nextpow2(L);
    %convert FFT bins to frequences
    f = 0:(Fs/N):(Fs/2);
end
