function Data = remove_outliers(Data)    
    maxval=0.25; % = 25mV; %recordings should be within [-25mV,25mV] after substracting baseling
    ids=find(abs(Data(:,1))>maxval);
    Data(ids,1)=sign(Data(ids,1))*maxval;
end