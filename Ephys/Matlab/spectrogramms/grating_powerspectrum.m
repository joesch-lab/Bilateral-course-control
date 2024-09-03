function grating_powerspectrum(grating_resfile)
% grating_resfile='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpD\HS\HSN\210623\res\stimuli_grating_stimuli_2021_06_23_17_27_43_raw_ds.mat';
% clear all;
% grating_resfile='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpND\HS\HSN\201202\res\stimuli_grating_stimuli_2020_12_02_16_37_24_raw_ds.mat';
load(grating_resfile);
sampling_rate = 10000;
if contains(exp_data.cell_grp, 'HS')
    pd_dir=1;
    nd_dir=5;
else
    pd_dir=3;
    nd_dir=7;
end

pd_trace=squeeze(exp_data.av_responses(1,pd_dir,:));
bl1_trace=squeeze(exp_data.av_baseline(1,pd_dir,:));
nd_trace=squeeze(exp_data.av_responses(1,nd_dir,:));
bl2_trace=squeeze(exp_data.av_baseline(1,nd_dir,:));
bl1mean=mean(bl1_trace);
bl2mean=mean(bl2_trace);

%normalize to mean baseline at zero;
bl1_trace = bl1_trace-bl1mean;
bl2_trace = bl2_trace-bl2mean;
combined_baseline = [bl1_trace; bl2_trace];

pd_trace = pd_trace-bl1mean;
nd_trace = nd_trace-bl2mean;

total_trace=[bl1_trace; pd_trace; bl2_trace; nd_trace];

%variance of indiv traces
bl1sum=sum(bl1_trace.^2);
bl2sum=sum(bl2_trace.^2);
pdsum=sum(pd_trace.^2);
ndsum=sum(nd_trace.^2);
combined_baseline_sum = sum(combined_baseline.^2);
totalsum = sum(total_trace.^2);

%find spectrum of total trace
tlen=length(total_trace);
t = (1:tlen)/sampling_rate;
xTable = timetable(seconds(t)', total_trace);    
[ptotal,ftotal] = pspectrum(xTable);
ptotalsum=sum(ptotal.^2);
ptotal_normed=ptotal/(ptotalsum.^0.5);

%find the index of the frequency=200
fcut=201;
fcut_ind = find(ftotal<fcut,1,'last');

%find spectrum of baseline 1
bl1len=length(bl1_trace);
t = (1:bl1len)/sampling_rate;
xTable = timetable(seconds(t)', bl1_trace);    
[pbl1,fbl1] = pspectrum(xTable);
pbl1_normed=(pbl1/(sum(pbl1.^2).^0.5))*(bl1sum/totalsum);
fcut_ind_bl1 = find(fbl1<fcut,1,'last');

%find spectrum of baseline 2
bl2len=length(bl2_trace);
t = (1:bl2len)/sampling_rate;
xTable = timetable(seconds(t)', bl2_trace);      
[pbl2,fbl2] = pspectrum(xTable);
pbl2_normed=(pbl2/(sum(pbl2.^2).^0.5))*(bl2sum/totalsum);
fcut_ind_bl2 = find(fbl2<fcut,1,'last');


%find spectrum of PD trace
pdlen=length(pd_trace);
t = (1:pdlen)/sampling_rate;
xTable = timetable(seconds(t)', pd_trace);   
[ppd,fpd] = pspectrum(xTable);
ppd_normed=(ppd/(sum(ppd.^2).^0.5))*(pdsum/totalsum);
fcut_ind_pd = find(fpd<fcut,1,'last');

%find spectrum of ND trace
ndlen=length(nd_trace);
t = (1:ndlen)/sampling_rate;
xTable = timetable(seconds(t)', nd_trace);   
[pnd,fnd] = pspectrum(xTable);
pnd_normed=(pnd/(sum(pnd.^2).^0.5))*(ndsum/totalsum);
fcut_ind_nd = find(fnd<fcut,1,'last');

%find spectrum of combined baseline
ndlen=length(combined_baseline);
t = (1:ndlen)/sampling_rate;
xTable = timetable(seconds(t)', combined_baseline);   
[pbl,fbl] = pspectrum(xTable);
pbl_normed=(pbl/(sum(pbl.^2).^0.5))*(combined_baseline_sum/totalsum);
fcut_ind_bl = find(fnd<fcut,1,'last');


%save all data
power_data.resfile = grating_resfile;
power_data.cell_grp = exp_data.cell_grp;
power_data.cell_type = exp_data.cell_type;
power_data.strain_type = exp_data.strain_type;
power_data.record_name = [exp_data.strain_type, ' ', exp_data.cell_type, ' ', datestr(exp_data.timestamp,'yyyy-mm-dd')];

power_data.traces.pd_trace=pd_trace;
power_data.traces.nd_trace=pd_trace;
power_data.traces.bl1_trace=bl1_trace;
power_data.traces.bl2_trace=bl2_trace;
power_data.traces.total_trace=total_trace;
power_data.traces.baseline=combined_baseline;
power_data.traces.meanval_bl1 = bl1mean;
power_data.traces.meanval_bl2 = bl2mean;
power_data.sumsignal.bl1 = bl1sum;
power_data.sumsignal.bl2 = bl2sum;
power_data.sumsignal.pd = pdsum;
power_data.sumsignal.nd = ndsum;
power_data.sumsignal.baseline = combined_baseline_sum;
power_data.sumsignal.total = totalsum;

power_data.powerspectrum.total.p= ptotal;
power_data.powerspectrum.total.f= ftotal;
power_data.powerspectrum.total.p_normed=ptotal_normed;
power_data.powerspectrum.total.f_cut_ind = fcut_ind;

power_data.powerspectrum.pd.p = ppd;
power_data.powerspectrum.pd.f= fpd;
power_data.powerspectrum.pd.p_normed=ppd_normed;
power_data.powerspectrum.pd.f_cut_ind = fcut_ind_pd;

power_data.powerspectrum.pd.pbl = pbl1;
power_data.powerspectrum.pd.fbl= fbl1;
power_data.powerspectrum.pd.pbl_normed=pbl1_normed;
power_data.powerspectrum.pd.fbl_cut_ind = fcut_ind_bl1;

power_data.powerspectrum.nd.p = pnd;
power_data.powerspectrum.nd.f= fnd;
power_data.powerspectrum.nd.p_normed=pnd_normed;
power_data.powerspectrum.nd.f_cut_ind = fcut_ind_nd;

power_data.powerspectrum.nd.pbl = pbl2;
power_data.powerspectrum.nd.fbl= fbl2;
power_data.powerspectrum.nd.pbl_normed=pbl2_normed;
power_data.powerspectrum.nd.fbl_cut_ind = fcut_ind_bl2;

power_data.powerspectrum.baseline.p = pbl;
power_data.powerspectrum.baseline.f= fbl;
power_data.powerspectrum.baseline.p_normed=pbl_normed;
power_data.powerspectrum.baseline.f_cut_ind = fcut_ind_bl;


spectrum_file = [grating_resfile(1:end-10),'pwspectrum.mat'];
save(spectrum_file, 'power_data');

% figure, 
% ind = power_data.powerspectrum.total.f_cut_ind;
% plot(power_data.powerspectrum.total.f(1:ind),power_data.powerspectrum.total.p_normed(1:ind));
% ind = power_data.powerspectrum.pd.f_cut_ind;
% hold on, plot(power_data.powerspectrum.pd.f(1:ind),power_data.powerspectrum.pd.p_normed(1:ind));
% ind = power_data.powerspectrum.nd.f_cut_ind;
% hold on, plot(power_data.powerspectrum.nd.f(1:ind),power_data.powerspectrum.nd.p_normed(1:ind));
% ind = power_data.powerspectrum.baseline.f_cut_ind;
% hold on, plot(power_data.powerspectrum.baseline.f(1:ind),power_data.powerspectrum.baseline.p_normed(1:ind));
% legend({"Total trace", "PD", "ND", "baseline"});
% xlabel('frequencies');
% ylabel('power');
% xlim([0,100]);
% title(power_data.record_name);
end







