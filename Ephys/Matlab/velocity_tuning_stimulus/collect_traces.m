% %collect velocity tuning curves, make plots per fly and average per strain
% clear all;
% celltype='HSN';
% strain1='FlpD';
% strain2='FlpND';
% 
% save_folder = 'C:\DATA\EPhys\vel_tuning';
% repo = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\';
% 
% %collect vel tuning log files for the strain and cell type
% log1set_name = fullfile(repo,'**\stimuli_grating_stimuli_HS_Vel_thick*.log');
% log1set=dir(log1set_name);
% %keep only file containing the celltype and the strain
% nf=length(log1set);
% todel=[];
% for i=1:nf
%     if ~(contains(log1set(i).folder,celltype) && contains(log1set(i).folder,strain1))
%         todel=[todel,i];
%     end    
% end
% log1set(todel)=[];
% 
% %collect vel tuning log files for the strain and cell type
% log2set_name = fullfile(repo,'**\stimuli_grating_stimuli_HS_Vel_thick*.log');
% log2set=dir(log2set_name);
% %keep only file containing the celltype and the strain
% nf=length(log2set);
% todel=[];
% for i=1:nf
%     if ~(contains(log2set(i).folder,celltype) && contains(log2set(i).folder,strain2))
%         todel=[todel,i];
%     end    
% end
% log2set(todel)=[];
% 
% 
% fig = figure;
% 
% speeddir_all=[0.1, 0.2, 0.5, 1, 2.5, 5, 7.5, 15];
% nspeeds_all=length(speeddir_all);
% sample_rate=10000;
% traces1_pd=[];%nan(1,nspeeds_all,3*sample_rate);
% traces1_nd=[];%nan(1,nspeeds_all,3*sample_rate);
% traces2_pd=[];%nan(1,nspeeds_all,3*sample_rate);
% traces2_nd=[];%nan(1,nspeeds_all,3*sample_rate);
% 
% 
% %process files with and collect the traces
% nf=length(log1set);
% for i=1:nf
%     fulllogfile=fullfile(log1set(i).folder, log1set(i).name);
%     %check if the file was analyzed
%     matfile = fullfile(log1set(i).folder,'res');
%     matfile_name = fullfile(matfile,'vel_tuning_*spf01.mat');
%     matinfo=dir(matfile_name);
%     disp([num2str(i),': Analyzing folder ',log1set(i).folder]);
%     save_velocity_tuning_raw_data_folder(log1set(i).folder);     
%     matfile=fullfile(matinfo(1).folder, matinfo(1).name);
%     load(matfile);
%     
%     mean_traces=squeeze(mean(exp_data.raw_traces,2));
%     speed_vec=exp_data.speed_dir(1:2:end,1)';
%     nspeeds=length(speed_vec); 
%     cols=brewermap(nspeeds,'RdYlBu');
%     cols=flipud(cols);
% 
%     %make figure with raw traces per file    
%     clf(fig);
%     ax1 = subplot(2,1,1);
%     x=1:length(exp_data.raw_traces);
%     x=x/sample_rate;
%     for j=1:nspeeds
%         plot(x, mean_traces((j-1)*2+1,:),'Color',cols(j,:));
%         hold on;
%     end    
%     legend(num2str(speed_vec'));
%     
%     ax2 = subplot(2,1,2);
%     for j=1:nspeeds
%         plot(x, mean_traces(j*2,:),'Color',cols(j,:));
%         hold on;
%     end
%     linkaxes([ax1,ax2],'x');
%     %take the last part of the folder path
%     [subpth,lastsubfolder,~]=fileparts(log1set(i).folder);
%     %if the last subfolder are all digits, it is date
% 
%     titlestr = [strain1, ' ', celltype, ' '];
%     exp='\d{6}';
%     if regexp(lastsubfolder,exp)
%         date=lastsubfolder;
%         titlestr =[titlestr, date];
%     else
%         animal_id=lastsubfolder;
%         [~,date,~]=fileparts(subpth);
%         titlestr =[titlestr, date,' ', animal_id];
%     end
%     suptitle(titlestr);
%     filename= [strrep(titlestr,' ','_'),'.fig'];
%     filename=fullfile(save_folder,filename);
%     savefig(fig,filename);
%     
%     
%     %normalize the traces by max and save into the average trace   
%     maxval=max(abs(mean_traces(:)));
%     nsamples=min(30000, size(mean_traces,2));
%       
%     speeds_here=nan(nspeeds,1);
%     for j=1:nspeeds
%         speeds_here(j)=find(speeddir_all==speed_vec(j));
%     end
%     pd_traces=mean_traces(1:2:end,:);
%     nd_traces=mean_traces(2:2:end,:);
%     this_traces=nan(nspeeds_all,30000);
%     this_traces(speeds_here,1:nsamples)=pd_traces/maxval;
%     traces1_pd=cat(3,traces1_pd,this_traces);    
%     this_traces(speeds_here,1:nsamples)=nd_traces/maxval;
%     traces1_nd=cat(3,traces1_nd,this_traces);    
% end
% 
% %repeat for the strain 2
% nf=length(log2set);
% for i=1:nf
%     fulllogfile=fullfile(log2set(i).folder, log2set(i).name);
%     %check if the file was analyzed
%     matfile = fullfile(log2set(i).folder,'res');
%     matfile_name = fullfile(matfile,'vel_tuning_*spf01.mat');
%     matinfo=dir(matfile_name);
%     disp([num2str(i),': Analyzing folder ',log2set(i).folder]);
%     save_velocity_tuning_raw_data_folder(log2set(i).folder);     
%     matfile=fullfile(matinfo(1).folder, matinfo(1).name);
%     load(matfile);
%         
%     mean_traces=squeeze(mean(exp_data.raw_traces,2));
%     speed_vec=exp_data.speed_dir(1:2:end,1)';
%     nspeeds=length(speed_vec); 
%     cols=brewermap(nspeeds,'RdYlBu');
%     cols=flipud(cols);
% 
%     %make figure with raw traces per file    
%     clf(fig);
%     ax1 = subplot(2,1,1);
%     x=1:length(exp_data.raw_traces);
%     x=x/sample_rate;
%     for j=1:nspeeds
%         plot(x, mean_traces((j-1)*2+1,:),'Color',cols(j,:));
%         hold on;
%     end    
%     legend(num2str(speed_vec'));
%     
%     ax2 = subplot(2,1,2);
%     for j=1:nspeeds
%         plot(x, mean_traces(j*2,:),'Color',cols(j,:));
%         hold on;
%     end
%     linkaxes([ax1,ax2],'x');
%     %take the last part of the folder path
%     [subpth,lastsubfolder,~]=fileparts(log2set(i).folder);
%     %if the last subfolder are all digits, it is date
% 
%     titlestr = [strain2, ' ', celltype, ' '];
%     exp='\d{6}';
%     if regexp(lastsubfolder,exp)
%         date=lastsubfolder;
%         titlestr =[titlestr, date];
%     else
%         animal_id=lastsubfolder;
%         [~,date,~]=fileparts(subpth);
%         titlestr =[titlestr, date,' ', animal_id];
%     end
%     suptitle(titlestr);
%     filename= [strrep(titlestr,' ','_'),'.fig'];
%     filename=fullfile(save_folder,filename);
%     savefig(fig,filename);
%     
%     
%     %normalize the traces by max and save into the average trace   
%     maxval=max(abs(mean_traces(:)));
%     nsamples=min(30000, size(mean_traces,2));
%       
%     speeds_here=nan(nspeeds,1);
%     for j=1:nspeeds
%         speeds_here(j)=find(speeddir_all==speed_vec(j));
%     end
%     pd_traces=mean_traces(1:2:end,:);
%     nd_traces=mean_traces(2:2:end,:);
%     this_traces=nan(nspeeds_all,30000);
%     this_traces(speeds_here,1:nsamples)=pd_traces/maxval;
%     traces2_pd=cat(3,traces2_pd,this_traces);    
%     this_traces(speeds_here,1:nsamples)=nd_traces/maxval;
%     traces2_nd=cat(3,traces2_nd,this_traces);    
% end
% 
% 
% traces1_pd=mean(traces1_pd,3);
% traces1_nd=mean(traces1_nd,3);
% traces2_pd=mean(traces2_pd,3);
% traces2_nd=mean(traces2_nd,3);
%     
    

traces1_pd_n=mean(traces1_pd,3,'omitnan');
traces1_nd_n=mean(traces1_nd,3,'omitnan');
traces2_pd_n=mean(traces2_pd,3,'omitnan');
traces2_nd_n=mean(traces2_nd,3,'omitnan');

figure,

subplot(3,4,[1,2])
plot(traces1_pd_n(1,1:10000),'--','Color', cols(end,:));
hold on;
plot(traces2_pd_n(1,1:10000),'Color', cols(end,:));
hold on;
plot(traces1_nd_n(1,1:10000),'--','Color', cols(1,:));
hold on;
plot(traces2_nd_n(1,1:10000),'Color', cols(1,:));
title(['Speed ',num2str(speeddir_all(1))]);


subplot(3,4,[5,6])
plot(traces1_pd_n(2,1:10000),'--','Color', cols(end,:));
hold on;
plot(traces2_pd_n(2,1:10000),'Color', cols(end,:));
hold on;
plot(traces1_nd_n(2,1:10000),'--','Color', cols(1,:));
hold on;
plot(traces2_nd_n(2,1:10000),'Color', cols(1,:));
title(['Speed ',num2str(speeddir_all(2))]);

subplot(3,4,[9,10])
plot(traces1_pd_n(3,1:10000),'--','Color', cols(end,:));
hold on;
plot(traces2_pd_n(3,1:10000),'Color', cols(end,:));
hold on;
plot(traces1_nd_n(3,1:10000),'--','Color', cols(1,:));
hold on;
plot(traces2_nd_n(3,1:10000),'Color', cols(1,:));
title(['Speed ',num2str(speeddir_all(3))]);


subplot(3,4,[3,4])
plot(traces1_pd_n(end-2,1:10000),'--','Color', cols(end,:));
hold on;
plot(traces2_pd_n(end-2,1:10000),'Color', cols(end,:));
hold on;
plot(traces1_nd_n(end-2,1:10000),'--','Color', cols(1,:));
hold on;
plot(traces2_nd_n(end-2,1:10000),'Color', cols(1,:));
title(['Speed ',num2str(speeddir_all(end-2))]);


subplot(3,4,[7,8])
plot(traces1_pd_n(end-1,1:10000),'--','Color', cols(end,:));
hold on;
plot(traces2_pd_n(end-1,1:10000),'Color', cols(end,:));
hold on;
plot(traces1_nd_n(end-1,1:10000),'--','Color', cols(1,:));
hold on;
plot(traces2_nd_n(end-1,1:10000),'Color', cols(1,:));
title(['Speed ',num2str(speeddir_all(end-1))]);

subplot(3,4,[11,12])
plot(traces1_pd_n(end,1:10000),'--','Color', cols(end,:));
hold on;
plot(traces2_pd_n(end,1:10000),'Color', cols(end,:));
hold on;
plot(traces1_nd_n(end,1:10000),'--','Color', cols(1,:));
hold on;
plot(traces2_nd_n(end,1:10000),'Color', cols(1,:));
title(['Speed ',num2str(speeddir_all(end))]);

suptitle('Velocity tuning traces: -- FlpD     - FlpND');




