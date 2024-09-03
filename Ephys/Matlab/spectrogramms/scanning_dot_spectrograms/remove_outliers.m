function Data = remove_outliers(Data)    
    maxval=0.2; % = 20mV; %recordings should be within [-20mV,20mV] after substracting baseling
    ids=find(abs(Data(:,1))>maxval);
    Data(ids,1)=sign(Data(ids,1))*maxval;
end