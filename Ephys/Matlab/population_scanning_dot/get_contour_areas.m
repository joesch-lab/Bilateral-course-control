function [level_area, level_cumarea] = get_contour_areas(cm)
%cm is the countour matrix returned by the matlab function 
%the function goes thru the list of of contours and computes the areas

% list of pairs: level - area, if there are several components of the same
% level, there area of each of them will be computed separately
level_area = [];

i=1;
while i<=size(cm,2)
leveli=cm(1,i);
npts=cm(2,i);
i=i+1;
com_cmi = mean(cm(:,i:i+npts-1),2)';
areai=0;
for pi=0:npts-2
    p1=cm(:,i+pi)';
    p2=cm(:,i+pi+1)';
    %area of the triangle
    ai=0.5*det([p1,1;p2,1;com_cmi,1]);
    areai = areai + abs(ai);
end
level_area=[level_area;[leveli,areai]];
i=i+npts;
end

%for each level compute cumulative area (sum area of all components)
ulevels=unique(level_area(:,1));
level_cumarea=zeros(length(ulevels),2);
for i=1:length(ulevels)
    uids=find(level_area(:,1)==ulevels(i));
    sum_a=sum(level_area(uids,2));
    level_cumarea(i,:)=[ulevels(i),sum_a];
end
end