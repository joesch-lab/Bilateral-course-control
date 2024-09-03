function [Data, Text_header, filename, samplingRate, minV, maxV]=openpr_flatten(filename, plot_yn, meanwin)
% this function open pr files form Labview and safes the data in a 2Dim
% array

if nargin < 1 || isempty(filename)
    error('myApp:argChk', 'No File Selected');
end


if ~exist('meanwin','var') || isempty(meanwin)
    meanwin=60;
end

fid = fopen(filename);

headervalue_01= fread(fid, [1], 'uint32','b');
values_02= fread(fid, [1,2], 'uint16','b');
AI_MAXMIN_03= fread(fid, [1,2], 'single','b');
SamplesRate_04= fread(fid, [1,3], 'uint32','b');
ChannelIndex_05 = fread(fid, [1], 'uint16','b');
tolong = 26 + 2*SamplesRate_04(3)+8*25;
Text_header = fread(fid, headervalue_01-tolong, '*char');

frewind(fid);
fread(fid, headervalue_01);
Datastring= fread(fid,'1000*double','b');
fclose(fid);
Data=zeros(length(Datastring)./2,2);

for i = 1: size(Datastring,1)./2000
    Data_temp(:,1) = Datastring((2*i-2)*1000+1:(2*i-1)*1000);
    Data_temp(:,2) = Datastring((2*i-2)*1000+1001:(2*i-1)*1000+1000);
    Data((i-1)*1000+1:(i)*1000,:)=Data_temp;
end

%find min and max values of the original Data
minV=min(Data(:,1));
maxV=max(Data(:,1));

%make sure the mean over large window is zero
samplingRate=SamplesRate_04(1);
meanwin=meanwin*samplingRate; %mean window
mm = movmean(Data(:,1),meanwin);
Data(:,1)=Data(:,1)-mm;

if plot_yn
    figure, plot(Data(:,2)./25,'b');
    hold on
    plot(Data(:,1),'k');
    hold off    
end
