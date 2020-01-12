
function [ cablecen ] = calculate_cable_centroid(cable,c42max,N)
%CALCULATE_CABLE_CENTROID Calculate centroid of points on a 2D domain with ...periodic boundary conditions
%   A = cable_centroid(array,c42cen). Where array has points ordered as [x1 y1 x2 y2 x3 y3 ..]
%   c42cen is the center of the GTP-Cdc42 patch

%store non-zero values in two separate matrices
actin_cable_X = [];
actin_cable_Y = [];

    for a=1:numel(cable(1,:))
        if (mod(a,2)== 1 && cable(1,a)~=0)
            actin_cable_X(round(a/2)) =  cable(1,a);
        elseif (mod(a,2)== 0 && cable(1,a)~=0)
            actin_cable_Y(a/2) = cable(1,a);
        end
    end
    
    %Find centroids by centering the patch
    
    %put everything near the center
    recenter_by = [50,50] - c42max(1,:) ;
    actin_cable_row_new = mod(actin_cable_X + recenter_by(1,1),N+10^-12) ; %mod(i-1,100)+1 to get values between 1 and 100
    actin_cable_col_new = mod(actin_cable_Y + recenter_by(1,2),N+10^-12);
    
    %find actin cable centroid
    cablecen = [mean(actin_cable_row_new) mean(actin_cable_col_new)]; %this should work as all cables are near the center now
    
    %move the center back to its original location
    cablecen = mod(cablecen - recenter_by - 1, 100) + 1;
    

end