function [spherical_coords, deformed_grid, deformedRF, ver_degree_span_half] = spherical_deformation(RF,hor_degree_span_half)
%the function stretches the image of RFx, RFy accross the sphere, such that
%the width of the images covers twice the hor_degree_span_half of the
%sphere

%the location of array elements in RF are assumed to be on the regular grid
% the deformation return the x and y coordinates of the deformed grid
deformedRF=zeros(size(RF));
deformed_grid=zeros(size(RF));
spherical_coords=zeros(size(RF));

[nrow, ncol,~] = size(RF);
ver_degree_span_half = 342*hor_degree_span_half/608;
azvals=linspace(-hor_degree_span_half,hor_degree_span_half,ncol);
elvals=linspace(-ver_degree_span_half,ver_degree_span_half,nrow);
xvals=1:ncol;
yvals=1:nrow;
[xmesh,ymesh]=meshgrid(xvals,yvals);
[azmesh,elmesh] = meshgrid(azvals,elvals);

for i=1:nrow
    for j=1:ncol
        azi=pi*azmesh(i,j)/180;
        eli=pi*elmesh(i,j)/180;
        [z,y3d,x3d]=sph2cart(eli,azi, 1);
        spherical_coords(i,j,:)=[round(azi*180/pi),round(eli*180/pi)];
        
        xn=xmesh(i,j)+RF(i,j,1);
        yn=ymesh(i,j)+RF(i,j,2);
        azi=pi*interp1(xvals,azvals,xn,"linear","extrap")/180;
        eli=pi*interp1(yvals,elvals,yn,"linear","extrap")/180;        
        [z,yn3d,xn3d]=sph2cart(eli,azi, 1);
        deformedRF(i,j,:)=[xn3d-x3d,yn3d-y3d];
        deformed_grid(i,j,:)=[x3d,y3d];
    end
end


end