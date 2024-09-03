function save_scanning_rect_powerspectrum_folder(curfolder)
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
                
        rawdata_i=trace_scanning_rect_one_recording_left_part(fulllog_file); 
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
    end
    
    Fs = 10000;          %sampling rate
    T = 1/Fs;           %sampling interval
    N = Fs/2*1000;%4096;%2048;%1024;      %FFT points
    nall=0;
    %convert FFT bins to frequences
    f = (0:N/2-1)*(Fs/N);
    f10=find(f<=10,1,'last');
    powerspectrum_av=zeros(1,numel(f));
    powerspectrum=zeros(total_nrep,numel(f));
    powfrac_low = zeros(total_nrep,1);
    powfrac_lowfft = zeros(total_nrep,1);
        
    %% for each rep and recording make a powerspectrum and average
    for i=1:ncount
        nrepi = size(data(i).trace,2);
        for ri=1:nrepi
            reptrace = data(i).trace(ri).reptrace;
            L = numel(reptrace); %Number of time points
            t = 0:T:(L-1)*T;    %time vector
            % FFT
            X = fft(reptrace,N);
            %make one sided signal
            SSB = X(1:N/2);
            SSB(2:end) = 2*SSB(2:end);            
            % Amplitude
            psp = abs(SSB/L); 
            
            powerspectrum_av=powerspectrum_av+psp;
            nall=nall+1;
            powerspectrum(nall,:)=psp;

            % spectrogramm to count HP regions
            [spow,fvec,tvec,ps] = spectrogram(reptrace,1024,512,0:1:250,10000,'yaxis');        
            sabs=abs(spow);
            powfrac_low(nall)= sum(sabs(1:11,:),'all')/sum(sabs,'all');
            powfrac_lowfft(nall)=sum(psp(1:f10))/sum(psp);
        end
    end

    powerspectrum_av=powerspectrum_av/nall;
                
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
    
    if isempty(logiles)
        return;
    end

    
    %create data structure    
    exp_data.folder=curfolder;
    exp_data.logfiles=logiles;
    exp_data.prfiles=prfiles;
    
    exp_data.pw.nrep= total_nrep;
    exp_data.pw.powerspectrum_all=powerspectrum;
    exp_data.pw.powerspectrum_av=powerspectrum_av;
    exp_data.pw.freq=f;
    exp_data.pw.powfrac_low=powfrac_low;
    exp_data.pw.powfrac_lowfft=powfrac_lowfft;
    
    exp_data.raw_trace=data;

    exp_data.cellinfo.strain = strain;
    exp_data.cellinfo.cellgroup = cellgroup;
    exp_data.cellinfo.celltype = celltype;
    exp_data.cellinfo.datestr = datestr;
    exp_data.cellinfo.cellid = cellid;       

    %save    
    if isempty(cellid)
        filename_i=['scaning_rect_pwspectrum_',strain,'_',celltype,'_',datestr];        
    else
        filename_i=['scaning_rect_pwspectrum_',strain,'_',celltype,'_',datestr,'_',cellid];        
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
