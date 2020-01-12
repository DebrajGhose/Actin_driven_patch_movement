% File name: endo_squeeze.m
% Squeezes membrane for exocytosis
% (exact) volume-preserving and (approx) shape-preserving interpolation
% Authors: Anita Layton and Daniel Lew
% Creation date: July 6, 2010
% fold: [input] function values at (xg,yg)
% fnew: [output] function values at (xold,yold)
% fwork: work space
% (patchx,patchy): indices where bin content set to 0
% scale: scaling factor needed to compute new volume

function [fnew] = endo_squeeze_radial ( fold, fwork, xg, yg, xold, yold, patchx, patchy, scale)
N = length(fold);

%approx shape-preserving interpolation
fwork(2:N+1,2:N+1) = fold;
fwork(1,2:N+1) = fold(N,1:N);
fwork(2:N+1,1) = fold(1:N,N);
fwork(1,1) = fold(N,N);

fnew = interp2(xg,yg,fwork,xold,yold);

%readjusting the central point to have a mean concentration 
x0 = [(N/2)-1:(N/2)+1];
y0 = [(N/2)-1:(N/2)+1];
% create regular grid
x0 = repmat(x0,3,1);
y0 = repmat(y0',1,3);

% x1 = [(N/2)-0.5:(N/2)+0.5];
% y1 = [(N/2)-0.5:(N/2)+0.5];
% % create regular grid
% x1 = repmat(x1,2,1);
% y1 = repmat(y1',1,2);

x1 = [(N/2)-0.5, (N/2)+0.5, N/2     ,  N/2     ];
y1 = [ N/2     ,  N/2     ,(N/2)-0.5, (N/2)+0.5];

fc = interp2(x0,y0,fold([(N/2)-1:(N/2)+1],[(N/2)-1:(N/2)+1]),x1,y1);
fnew(N/2,N/2) = mean(fc);

% adjust for volume conservation
total    = ((sum(sum(fold)) - sum(sum(fold(patchx,patchy)))))*scale;
new_tot  = sum(sum(fnew));
vol_diff = (total - new_tot);

% scaling to make sure new values don't exceed previous max and min
scaling  = (fold-min(min(fold))).*(max(max(fold))-fold);
fnew     = fnew + vol_diff.*scaling./sum(sum(scaling));

sc_tpm = vol_diff.*scaling./sum(sum(scaling));

if (~isempty(find(fnew<0)))
    
    neg_ind    = find(fnew < 0);            % finding negative grid points
    ind2       = setdiff([1:N^2],neg_ind);  % excluding them from mass 
                                            % redistribution
    mass_rd    = sum(fnew(neg_ind))/length(ind2);
    fnew(ind2) = fnew(ind2) + mass_rd;
    fnew(neg_ind) = 0;                      % setting negative bins to zero
    
end 


return

