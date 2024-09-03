function save_anystim_powerspectrum_folder(curfolder, stim_root)
% for the log file specified with 'stim_root' (eg stimuli_full_field_flashes) 
% in the cur folder, find the matching pr file and compute power spectrum analysis
% for the entire recording independetly on the parameters of the stimulus.

% this script will be used primary for analysing shakB2 recordings

% curfolder='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpND\VS\VS1_4\210128';
    %% get strain, cell type and date from the folder path
    [strain,cellgroup,celltype,datestr,cellid] = cell_info_from_path(curfolder);
    logpatt=['*',stim_root,'*.log'];
    logfiles=fullfile(curfolder,logpatt);
    logfileinfo=dir(logfiles);

    nfiles=length(logfileinfo);

    if nfiles==0
        return;
    end

    resfolder=fullfile(curfolder,'res');
    if ~exist(resfolder,'dir')
        mkdir(resfolder);
    end
    
    logiles={};
    prfiles={};

    Fs = 10000;          %sampling rate
    T = 1/Fs;           %sampling interval
    N = Fs/2*1000;%4096;%2048;%1024;      %FFT points
    nall=0;
    %convert FFT bins to frequences
    f = (0:N/2-1)*(Fs/N);
    f10=find(f<=10,1,'last');
    powerspectrum_av=zeros(1,numel(f));
    powerspectrum=zeros(nfiles,numel(f));
    powfrac_low = zeros(nfiles,1);
    powfrac_lowfft = zeros(nfiles,1);

    for fi=1:nfiles
        logfile_fullname=fullfile(logfileinfo(fi).folder,logfileinfo(fi).name);        
        prfilename = find_pr_file_closest_date(logfile_fullname);
        if isempty(prfilename)
            continue;
        end
    
        nall=nall+1;
        prfullname=fullfile(curfolder,prfilename);
        [Data, ~, ~, ~]=openpr_flatten_translate(prfullname,0);
        trace = Data(:,1)';        
            
        %structure to pool all data        
        logiles{nall}=logfile_fullname;        
        prfiles{nall}=prfilename;
        data(nall).trace=trace;      
            
        %% make a powerspectrum
        L = numel(trace); %Number of time points
        t = 0:T:(L-1)*T;    %time vector
        % FFT
        X = fft(trace,N);
        %make one sided signal
        SSB = X(1:N/2);
        SSB(2:end) = 2*SSB(2:end);            
        % Amplitude
        psp = abs(SSB/L);     
        
        powerspectrum_av=powerspectrum_av+psp;
        powerspectrum(nall,:)=psp;
    
        % spectrogramm to count HP regions
        [spow,fvec,tvec,ps] = spectrogram(trace,1024,512,0:1:250,10000,'yaxis');        
        sabs=abs(spow);
        powfrac_low(nall)= sum(sabs(1:11,:),'all')/sum(sabs,'all');
        powfrac_lowfft(nall)=sum(psp(1:f10))/sum(psp);
    end 
    powerspectrum_av=powerspectrum_av/nall;
    if nall<nfiles
        powerspectrum(nall+1:end,:)=[];
        powfrac_low(nall+1:end)=[];
        powfrac_lowfft(nall+1:end)=[];
    end
                
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
    
    exp_data.pw.nrep= nall;
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
        filename_i=['anystim_pwspectrum_',strain,'_',celltype,'_',datestr];        
    else
        filename_i=['anystim_pwspectrum_',strain,'_',celltype,'_',datestr,'_',cellid];        
    end    
    fullfile_i=fullfile(resfolder,[filename_i,'.mat']);
    save(fullfile_i,'exp_data','-v7.3');  
end



 %function to find the pr file closes to the log file     
function prfilename = find_pr_file_closest_date(logfilename)
    [folder,~,~]=fileparts(logfilename);    
    %date info from log file
    log_info = dir(logfilename);
    log_datevec=datevec(log_info.date);
    D = dir(folder);
    min_et=3600;
    prfilename=-1;
    %find the closest pr file by time
    for k = 3:length(D) % avoid using the first ones
       currfile = D(k).name; %file name
       if contains(currfile,'.pr','IgnoreCase', true)
           et=abs(etime(datevec(D(k).date),log_datevec));
           if et<min_et
               min_et=et;
               prfilename=currfile;
           end
       end
    end
    if min_et>10
       disp('No matching pr file has been found');
       prfilename=-1;
    end      
end
