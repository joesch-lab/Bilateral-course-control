function save_scanning_rect_powerspectrum_all_sides_folder(curfolder)
% for all scanning rect stimuli get the part of the trace
% corresponding to the time when rect was on the left
% for each such trace make a power spectrum and average it

% curfolder='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpND\VS\VS1_4\210128';
    %% get strain, cell type and date from the folder path
    [strain,cellgroup,celltype,datestr,cellid] = cell_info_from_path(curfolder);
    
    %% get all stimulus files
    filelist = dir(curfolder);%get list of files
    filelist = filelist(~[filelist.isdir]);  %remove folders from list

    stim_name='_scanning_rect_';
    Fs = 10000; %sampling rate

    %approximate duration of the traces: max 100s
    L_appr=100*Fs;
    % frequency bins
    f =  fft_bins(Fs,L_appr);
    % we don't care about frequencies higher than 1K
    fmax=10^3;
    fmaxind=find(f<=fmax,1,'last');
    f=f(1:fmaxind);

    
    %go thru the list and collect files with'_scanning_rect_' 
    stim_files={};
    i=1;
    for li =1:length(filelist)
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind) &&  contains(filelist(li).name,'.log','IgnoreCase', true)  
            stim_files{i}=filelist(li).name;            
            i=i+1;
        end
    end
    
    resfolder=fullfile(curfolder,'res');
    if ~exist(resfolder,'dir')
        mkdir(resfolder);
    end
    
    total_nrep=0;
    %% iterate thru stim files, get raw RF values
    nfiles=length(stim_files);
    ncount=0;
    logiles={};
    prfiles={};
    for i=1:nfiles       
       
        fulllog_file=fullfile(curfolder,stim_files{i});
        stimparam_i= read_scanning_rect_with_pause(fulllog_file);
        if stimparam_i.rheight~=30, continue; end
                
        rawdata_i=trace_scanning_rect_one_recording_left_right_gray(fulllog_file); 
        if isempty(rawdata_i)            
            continue;
        end
        ncount=ncount+1;
        %structure to pool all data        
        logiles{ncount}=stim_files{i};        
        prfiles{ncount}=rawdata_i.prname;
        nreps=stimparam_i.nrep;
        total_nrep=total_nrep+nreps;
        data(ncount).trace=rawdata_i.data.trace;        
        data(ncount).graytrace=rawdata_i.data.trace_gray;        
    end
    

    %% fft power spectrum for all traces
    powerspectrum_contra=zeros(total_nrep,fmaxind);    
    powerspectrum_contra_h=zeros(total_nrep,fmaxind);    
    powerspectrum_contra_v=zeros(total_nrep,fmaxind);
    powerspectrum_ipsi=zeros(total_nrep,fmaxind);
    powerspectrum_ipsi_h=zeros(total_nrep,fmaxind);
    powerspectrum_ipsi_v=zeros(total_nrep,fmaxind);
    powerspectrum_gray = zeros(nfiles,fmaxind);
    
    %contribution of parts of spectrum: low, mid, total
    pow_conra = zeros(total_nrep,3);    
    pow_conra_h = zeros(total_nrep,3);
    pow_conra_v = zeros(total_nrep,3);
    pow_ipsi = zeros(total_nrep,3);
    pow_ipsi_h= zeros(total_nrep,3);
    pow_ipsi_v = zeros(total_nrep,3);
    pow_gray = zeros(nfiles,3);
    
    nall=1;    
    for i=1:ncount
        nrepi = size(data(i).trace,2);
        for re=1:nrepi
            temptr=data(i).trace.reptrace_contra;
            [powerspectrum_contra(nall,:), pow_conra(nall,1), pow_conra(nall,2), pow_conra(nall,3),~] = get_trace_pwspectrum_maxlen(temptr, Fs, L_appr, fmax);
            temptr=data(i).trace.reptrace_contra_h;
            [powerspectrum_contra_h(nall,:), pow_conra_h(nall,1), pow_conra_h(nall,2), pow_conra_h(nall,3),~] = get_trace_pwspectrum_maxlen(temptr, Fs, L_appr, fmax);
            temptr=data(i).trace.reptrace_contra_v;
            [powerspectrum_contra_v(nall,:), pow_conra_v(nall,1), pow_conra_v(nall,2), pow_conra_v(nall,3),~] = get_trace_pwspectrum_maxlen(temptr, Fs, L_appr, fmax);
            temptr=data(i).trace.reptrace_ipsi;
            [powerspectrum_ipsi(nall,:), pow_ipsi(nall,1), pow_ipsi(nall,2), pow_ipsi(nall,3),~] = get_trace_pwspectrum_maxlen(temptr, Fs, L_appr, fmax);
            temptr=data(i).trace.reptrace_ipsi_h;
            [powerspectrum_ipsi_h(nall,:), pow_ipsi_h(nall,1), pow_ipsi_h(nall,2), pow_ipsi_h(nall,3),~] = get_trace_pwspectrum_maxlen(temptr, Fs, L_appr, fmax);
            temptr=data(i).trace.reptrace_ipsi_v;
            [powerspectrum_ipsi_v(nall,:), pow_ipsi_v(nall,1), pow_ipsi_v(nall,2), pow_ipsi_v(nall,3),~] = get_trace_pwspectrum_maxlen(temptr, Fs, L_appr, fmax);
            nall=nall+1; 
        end
        [powerspectrum_gray(i,:), pow_gray(nall,1), pow_gray(nall,2), pow_gray(nall,3),~] = get_trace_pwspectrum_maxlen(data(i).graytrace, Fs, L_appr, fmax);        
    end
    
    %make averages of the powespectra
    powerspectrum_contra_av=mean(powerspectrum_contra,1);    
    powerspectrum_contra_h_av=mean(powerspectrum_contra_h,1);        
    powerspectrum_contra_v_av=mean(powerspectrum_contra_v,1);    
    powerspectrum_ipsi_av=mean(powerspectrum_ipsi,1);    
    powerspectrum_ipsi_h_av=mean(powerspectrum_ipsi_h,1);
    powerspectrum_ipsi_v_av=mean(powerspectrum_ipsi_v,1);
    powerspectrum_gray_av=mean(powerspectrum_gray,1);
    
    %contribution of parts of spectrum: low, mid, total
    pow_conra_av = mean(pow_conra,1);  
    pow_conra_h_av = mean(pow_conra_h,1);
    pow_conra_v_av = mean(pow_conra_v,1);
    pow_ipsi_av = mean(pow_ipsi,1);
    pow_ipsi_h_av= mean(pow_ipsi_h,1);
    pow_ipsi_v_av = mean(pow_ipsi_v,1);
    pow_gray_av = mean(pow_gray,1);

    figure, 
    loglog(f, powerspectrum_contra_v_av); hold on; % no stim
    loglog(f, powerspectrum_ipsi_h_av); hold on; % stim
    loglog(f, powerspectrum_gray_av); hold on; % gray
    legend(["contra v","ipsi h","gray"]);
    title(strjoin({strain,cellgroup,celltype,datestr,cellid},' '),'Interpreter','none');
    
    
    %create data structure    
    exp_data.folder=curfolder;
    exp_data.logfiles=logiles;
    exp_data.prfiles=prfiles;
    
    exp_data.pw.nrep= total_nrep;   
    exp_data.pw.nfiles = nfiles;   
    exp_data.pw.freq=f;   

    %% save powespectra from all reps
    exp_data.pw_raw.powerspectrum_contra=powerspectrum_contra;
    exp_data.pw_raw.powerspectrum_contra_h=powerspectrum_contra_h;
    exp_data.pw_raw.powerspectrum_contra_v=powerspectrum_contra_v;
    exp_data.pw_raw.powerspectrum_ipsi=powerspectrum_ipsi;
    exp_data.pw_raw.powerspectrum_ipsi_h=powerspectrum_ipsi_h;
    exp_data.pw_raw.powerspectrum_ipsi_v=powerspectrum_ipsi_v;
    exp_data.pw_raw.powerspectrum_gray=powerspectrum_gray;

    exp_data.pw_raw.pow_conra=pow_conra;
    exp_data.pw_raw.pow_conra=pow_conra_h;
    exp_data.pw_raw.pow_conra_v=pow_conra_v;
    exp_data.pw_raw.pow_ipsi=pow_ipsi;
    exp_data.pw_raw.pow_ipsi_h=pow_ipsi_h;
    exp_data.pw_raw.pow_ipsi_v=pow_ipsi_v;
    exp_data.pw_raw.pow_gray=pow_gray;

    %% save average powespectra 
    exp_data.pw_av.powerspectrum_contra=powerspectrum_contra_av;
    exp_data.pw_av.powerspectrum_contra_h=powerspectrum_contra_h_av;
    exp_data.pw_av.powerspectrum_contra_v=powerspectrum_contra_v_av;
    exp_data.pw_av.powerspectrum_ipsi=powerspectrum_ipsi_av;
    exp_data.pw_av.powerspectrum_ipsi_h=powerspectrum_ipsi_h_av;
    exp_data.pw_av.powerspectrum_ipsi_v=powerspectrum_ipsi_v_av;
    exp_data.pw_av.powerspectrum_gray=powerspectrum_gray_av;

    exp_data.pw_av.pow_conra=pow_conra_av;
    exp_data.pw_av.pow_conra_h=pow_conra_h_av;
    exp_data.pw_av.pow_conra_v=pow_conra_v_av;
    exp_data.pw_av.pow_ipsi=pow_ipsi_av;
    exp_data.pw_av.pow_ipsi_h=pow_ipsi_h_av;
    exp_data.pw_av.pow_ipsi_v=pow_ipsi_v_av;
    exp_data.pw_av.pow_gray=pow_gray_av;

    %% save normalized average powespectra: the area under the curve is 1
    exp_data.pw_av_norm.powerspectrum_contra=powerspectrum_contra_av/sum(powerspectrum_contra_av);
    exp_data.pw_av_norm.powerspectrum_contra_h=powerspectrum_contra_h_av/sum(powerspectrum_contra_h_av);
    exp_data.pw_av_norm.powerspectrum_contra_v=powerspectrum_contra_v_av/sum(powerspectrum_contra_v_av);
    exp_data.pw_av_norm.powerspectrum_ipsi=powerspectrum_ipsi_av/sum(powerspectrum_ipsi_av);
    exp_data.pw_av_norm.powerspectrum_ipsi_h=powerspectrum_ipsi_h_av/sum(powerspectrum_ipsi_h_av);
    exp_data.pw_av_norm.powerspectrum_ipsi_v=powerspectrum_ipsi_v_av/sum(powerspectrum_ipsi_v_av);
    exp_data.pw_av_norm.powerspectrum_gray=powerspectrum_gray_av/sum(powerspectrum_gray_av);

    exp_data.pw_av_norm.pow_conra=pow_conra_av/pow_conra_av(3);
    exp_data.pw_av_norm.pow_conra_h=pow_conra_h_av/pow_conra_h_av(3);
    exp_data.pw_av_norm.pow_conra_v=pow_conra_v_av/pow_conra_v_av(3);
    exp_data.pw_av_norm.pow_ipsi=pow_ipsi_av/pow_ipsi_av(3);
    exp_data.pw_av_norm.pow_ipsi_h=pow_ipsi_h_av/pow_ipsi_h_av(3);
    exp_data.pw_av_norm.pow_ipsi_v=pow_ipsi_v_av/pow_ipsi_v_av(3);
    exp_data.pw_av_norm.pow_gray=pow_gray_av/pow_gray_av(3);
   
    %% save traces
    exp_data.raw_trace=data;
    %% cell info
    exp_data.cellinfo.strain = strain;
    exp_data.cellinfo.cellgroup = cellgroup;
    exp_data.cellinfo.celltype = celltype;
    exp_data.cellinfo.datestr = datestr;
    exp_data.cellinfo.cellid = cellid;       

    %save    
    if isempty(cellid)
        filename_i=['scaning_rect_powspectrum_',strain,'_',celltype,'_',datestr];        
    else
        filename_i=['scaning_rect_powspectrum_',strain,'_',celltype,'_',datestr,'_',cellid];        
    end    
    fullfile_i=fullfile(resfolder,[filename_i,'.mat']);
    save(fullfile_i,'exp_data','-v7.3');  
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

function f = fft_bins(Fs, L)  
    N=2^nextpow2(L);
%     N = Fs/2*1000;%2^nextpow2(L);
    %convert FFT bins to frequences
    f = 0:(Fs/N):(Fs/2);
end