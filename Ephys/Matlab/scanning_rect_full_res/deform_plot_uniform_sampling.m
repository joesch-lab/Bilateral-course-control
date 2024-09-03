function [bin_dx_av, bin_resp, bin_resp_deformed, bin_hor, bin_ver, bin_hor_sep, bin_ver_sep, el_bincenters, az_bincenters] = deform_plot_uniform_sampling(RFx,RFy, nbin_hor, resp_xyc, plot_yn)

    % parameters of the stimulus, usually the same for scanning dot and
    % brownian dots
    w=608;
    h=342;
    wall_offset = [90,20];
    wall_thickness=5;
    xoffset = wall_offset(1)+wall_thickness;
    yoffset = wall_offset(2)+wall_thickness;
    
    %assuming the stimulus convolved with neuronal activity is the full
    %size of the screen
%     h=size(RFx,1);
%     w=size(RFx,2);
    
    %number of bins horizontally
    nbinx=nbin_hor;
    %number of bins vertically defined by the same size of the bin
    nbiny=ceil((h-2*yoffset)/((w-2*xoffset)/nbinx));
%     nbiny=ceil(h/(w/nbinx));
      
    % %for debuging init some RF
    % RFx = ones(14,15)*0.7;
    % RFy = ones(14,15)*0.5;

    %size of the RF
    RF_y_size=size(RFx,1);
    RF_x_size=size(RFx,2);
    if resp_xyc=='x'
        RF_resp=RFx;
    elseif resp_xyc=='y'
        RF_resp=RFy;
    else
        RF_resp=sqrt(RFx.^2+RFy.^2);
    end
%     xoffset=0;
%     yoffset=0;
%     
    %how many pixels from the stimulus image the RF pixel spans    
    pix_size_x=(w-2*xoffset)/RF_x_size;
    pix_size_y=(h-2*yoffset)/RF_y_size;
%     pix_size_x=1;%(w-2*xoffset)/RF_x_size;
%     pix_size_y=1;%(h-2*yoffset)/RF_y_size;

    % x and y values of the pixels (centers)    
    x_or = arrayfun(@(i) xoffset+(i-1)*pix_size_x+0.5*pix_size_x, 1:RF_x_size);
    y_or = arrayfun(@(i) yoffset+(i-1)*pix_size_y+0.5*pix_size_y, 1:RF_y_size);
    %repeat columns and rows of the coordinates
    Xor = repelem(x_or,RF_y_size,1);
    Yor = repelem(y_or',1, RF_x_size);
    
    %coordinates of the end of the arrow in the RF
    Xnew = Xor + RFx;
    Ynew = Yor + RFy;

    %transform 2D indices to 3D coordinates
    RF3D=zeros(RF_y_size,RF_x_size,3); %3D location of grid points
    dx3D=zeros(RF_y_size,RF_x_size,3); %3D vectors of optic flow
    
    for i=1:RF_y_size
        for j=1:RF_x_size
            [x,y,z, el, az] = get_3D_coordinates(Xor(i,j),Yor(i,j));
            [xn,yn,zn, eln, azn] = get_3D_coordinates(Xnew(i,j),Ynew(i,j));            
            dx3D(i,j,:)=[xn-x,yn-y,zn-z];
            RF3D(i,j,:)=[x,y,z];           
        end
    end
    
    %plot only X and Y positions
    X_proj=RF3D(:,:,1); X_proj=X_proj(:);
    Y_proj=RF3D(:,:,2); Y_proj=Y_proj(:);
    Z_proj=RF3D(:,:,3); Z_proj=Z_proj(:);
    RF_resp=RF_resp(:);
%     figure, scatter(X_proj,Y_proj,'.');    
    
    Xdx=dx3D(:,:,1); Xdx=Xdx(:);
    Ydx=dx3D(:,:,2); Ydx=Ydx(:);
    Zdx=dx3D(:,:,3); Zdx=Zdx(:);
   
    Xdeform = X_proj;
    Ydeform = Y_proj;
    npts=length(Xdeform);
    Xrange=max(Xdeform)-min(Xdeform);
    Yrange=max(Ydeform)-min(Ydeform);
    minX=min(Xdeform);
    minY=min(Ydeform);    
   
    %make lower-res picture 
    %size of the lowres bin
    xbinsize=Xrange/nbinx;
    ybinsize=Yrange/nbiny;
    
    xbincenters=arrayfun(@(i) minX+xbinsize/2+(i-1)*xbinsize,1:nbinx);
    ybincenters=arrayfun(@(i) minY+ybinsize/2+(i-1)*ybinsize,1:nbiny);
    %azimuth and elevation values of bins in radians
    el_bincenters=ybincenters*pi/w;
    az_bincenters=xbincenters*pi/w;    
    
    %average the direction vector for each bin 
    bin_dx=zeros(nbiny,nbinx,3);
    bin_resp=zeros(nbiny,nbinx);    
    bin_count=zeros(nbiny,nbinx);
   
    %count and sum up high res points in the low res grid
    for i=1:npts
        %transform x coordinate to the bin id  
        bin_raw=(Xdeform(i)-minX)/xbinsize;          
        if floor(bin_raw)==0 %min values is in the 1st bin
            bin_xc=1;
        else
            bin_xc = min(ceil(bin_raw),nbinx);
        end
       
        %transform y coordinate to the bin id  
        bin_raw=(Ydeform(i)-minY)/ybinsize;            
        if floor(bin_raw)==0 %min values is in the 1st bin
            bin_yc=1;
        else
            bin_yc = min(ceil(bin_raw),nbiny);
        end
        
        bin_count(bin_yc,bin_xc)=bin_count(bin_yc,bin_xc)+1;
        bin_dx(bin_yc,bin_xc,:)=bin_dx(bin_yc,bin_xc,:)+reshape([Xdx(i),Ydx(i),Zdx(i)],1,1,3);       
        bin_resp(bin_yc,bin_xc)=bin_resp(bin_yc,bin_xc)+RF_resp(i);       
    end
    %sum up direction vectors horizontally
    bin_dx_hor=sum(bin_dx(:,:,1));
    bin_dy_hor=sum(bin_dx(:,:,2));
    bin_count_hor=sum(bin_count);
    bin_dx_hor=bin_dx_hor./bin_count_hor;
    bin_dx_hor(isnan(bin_dx_hor))=0;      
    bin_dy_hor=bin_dy_hor./bin_count_hor;
    bin_dy_hor(isnan(bin_dy_hor))=0;   
    bin_hor_sep=[bin_dx_hor;bin_dy_hor];
    
    %sum up direction vectors vertically
    bin_dx_ver=sum(bin_dx(:,:,1),2);
    bin_dy_ver=sum(bin_dx(:,:,2),2);
    bin_count_ver=sum(bin_count,2);
    bin_dx_ver=bin_dx_ver./bin_count_ver;
    bin_dx_ver(isnan(bin_dx_ver))=0;      
    bin_dy_ver=bin_dy_ver./bin_count_ver;
    bin_dy_ver(isnan(bin_dy_ver))=0;   
    bin_ver_sep=[bin_dx_ver';bin_dy_ver'];
    
    %make average on 2D grid
    bin_dx_av=bin_dx./bin_count;
    bin_dx_av(isnan(bin_dx_av))=0;  
    
    %length of RF vectors
    bin_resp_deformed = sqrt(bin_dx_av(:,:,1).^2+bin_dx_av(:,:,2).^2);
    bin_hor=mean(bin_resp_deformed,'omitnan');
    bin_ver=mean(bin_resp_deformed,2,'omitnan');
    bin_resp_deformed(isnan(bin_resp_deformed))=0;  

    %responses without taking into account deformations
    bin_resp=bin_resp./bin_count;
    bin_resp(isnan(bin_resp))=0;      

    % we care only about X and Y component of the RF:   
    bin_dx_av=bin_dx_av(:,:,1:2);
    
    %if the plot was requested
    if plot_yn
%         figure, bar(1:nbinx,[bin_dx_hor;bin_dy_hor]');
%         figure, barh(1:nbiny,[bin_dx_ver';bin_dy_ver']');
 
        half_az_spna_deg=70;
    %     figure, imagesc(bin_count);
        xi=1:nbinx;
        yi=1:nbiny;
        [Xplot,Yplot]= meshgrid(xi,yi);

        figure, quiver(Xplot,Yplot,bin_dx_av(:,:,1),bin_dx_av(:,:,2));    
        xlim([0,nbinx+1]);
        ylim([0,nbiny+1]);   
        set(gca,'XDir','reverse'); 
        xticks([1,nbinx/2,nbinx]);
        yticks([nbiny/2]);
        yticklabels(0);
        xlabels=[half_az_spna_deg,0,-half_az_spna_deg];
        xticklabels(xlabels);
        xlabel('azimuth');
        ylabel('elevation');
    end

% %     plot 3D quiver
%     X=RF3D(:,:,1); X=X(:);
%     Y=RF3D(:,:,2); Y=Y(:);
%     Z=RF3D(:,:,3); Z=Z(:);
%     dX=dx3D(:,:,1); dX=dX(:);
%     dY=dx3D(:,:,2); dY=dY(:);
%     dZ=dx3D(:,:,3); dZ=dZ(:);
% % 
%     figure, quiver3(X,Y,Z,dX, dY, dZ);
% %     figure, quiver3(Z,Y,X,dZ,dY,dX);
%     grid off
%     axis equal
%     view(0,90)
%     xlabel('x')
%     ylabel('y')
end



function [x,y,z, el, az] = get_3D_coordinates(xim,yim)

    h=342;
    w=608;
    % h2=h/2;
    % w2=w/2;
    % m2=max(h2,w2);
    m2=(w+1)/2;
    pi2=pi/2;
    % r=1;
    %radius of the sphere so that the long side of the image wraps half sphere without scaling up or down
    r=w/pi; 
    %padding for the shorter side of the image
    v_pad=floor((w-h)/2);

    el=(yim+v_pad-1-m2)*pi2/m2;
    az=(xim-m2)*pi2/m2;
    [z,y,x]=sph2cart(el,az, r); 
end


