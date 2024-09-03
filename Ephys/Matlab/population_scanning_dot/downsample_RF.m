function RF_lowres = downsample_RF(RFx, RFy, num_hor_bins, num_ver_bins)
    [nrow,ncol]=size(RFx);
    if ~exist('num_ver_bins','var') || isempty(num_ver_bins)
        num_ver_bins=ceil(nrow*num_hor_bins/ncol);
    end
    hor_edges=round(linspace(1,ncol,num_hor_bins+1));
    ver_edges=round(linspace(1,nrow,num_ver_bins+1));
    
    RF_lowres=zeros(num_ver_bins,num_ver_bins,2);
    %compute mean x,y responses and patch values
    for ri=1:num_ver_bins
        for ci=1:num_hor_bins
            subRF = RFx(ver_edges(ri):ver_edges(ri+1),hor_edges(ci):hor_edges(ci+1));
            RF_lowres(ri,ci,1)=mean(subRF(:));
            subRF = RFy(ver_edges(ri):ver_edges(ri+1),hor_edges(ci):hor_edges(ci+1));
            RF_lowres(ri,ci,2)=mean(subRF(:));
        end
    end
end