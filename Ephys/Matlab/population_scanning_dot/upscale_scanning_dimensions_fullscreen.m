 function [RFx_hres, RFy_hres]= upscale_scanning_dimensions_fullscreen(RFx,RFy)
    %find non-zero rows and columns
    sumrows=sum(RFx,2);
    sumcols=sum(RFx,1);
    nnzrows=find(sumrows~=0);
    nnzcols=find(sumcols~=0);
    RFxeffective=RFx(nnzrows,nnzcols);
    
    sumrows=sum(RFy,2);
    sumcols=sum(RFy,1);
    nnzrows=find(sumrows~=0);
    nnzcols=find(sumcols~=0);
    RFyeffective=RFy(nnzrows,nnzcols);    
    
%     %upscale non-zero components
%     maxrows=max(size(RFxeffective,1), size(RFyeffective,1));
%     maxcols=max(size(RFxeffective,2), size(RFyeffective,2));
    
    RFx_hres=imresize(RFxeffective,[342,608],'nearest');
    RFy_hres=imresize(RFyeffective,[342,608],'nearest');    
end
