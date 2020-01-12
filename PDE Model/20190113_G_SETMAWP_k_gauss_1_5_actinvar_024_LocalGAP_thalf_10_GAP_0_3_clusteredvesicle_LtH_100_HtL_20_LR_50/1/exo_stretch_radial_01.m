% fold: [input] function values at (xg,yg)          - regular mesh
% fnew: [output] function values at (xold,yold)     - transformed mesh
% fwork: work space
% (patchx,patchy): indices where bin content set to 0
% scale: scaling factor needed to compute new volume

function fnew = exo_stretch_radial_01 ( fold, fwork, xg, yg, xold, yold, patchx, patchy, scale)
N = length(fold);

% approx shape-preserving interpolation
fwork(2:N+1,2:N+1) = fold;
fwork(1,2:N+1) = fold(N,1:N);
fwork(2:N+1,1) = fold(1:N,N);
fwork(1,1) = fold(N,N);

fnew = interp2(xg,yg,fwork,xold,yold);

% set new patch of cells content to 0
fnew(patchx,patchy) = 0;

% adjust for volume conservation
total    = sum(sum(fold))*scale;        % old mass
new_tot  = sum(sum(fnew));              % new mass
vol_diff = total - new_tot;             % difference in mass

% don't distribute to new cells
                    % [indx,indy] = index of winning bin, bottom-left of 2x2 spot
                    % patchx      = [indx:indx+1]
                    % patchy      = [indy:indy+1]
my_patchx = repmat(patchx,1,length(patchx)); 
                    % my_patchx = [x1,x2,x1,x2]
my_patchy = reshape(repmat(patchy,length(patchx),1),1,length(patchx)^2);
                    % repmat(patchy,length(patchx),1) = [y1,y2;y1,y2]
                    % reshape([y1,y2;y1,y2],1,length(patchx)^2) = [y1,y1,y2,y2]
                    %
                    % my_patchx = [x1,x2,x1,x2]
                    % my_patchy = [y1,y1,y2,y2]                    
new_ind = sub2ind([N,N],my_patchx,my_patchy);
                    % sub2ind(matrixSize, rowSub, colSub)
                    % sub2ind([N,N],(x1,y1),(x2,y1),(x1,y2),(x2,y2))
ind = setdiff([1:N^2],new_ind);

% scaling to make sure new values don't exceed previous max and min
scaling = (fold-min(min(fold))).*(max(max(fold))-fold);
sc_tpm = vol_diff.*scaling./sum(sum(scaling));


% adding new mass as a concentration
fnew(ind) = fnew(ind) + vol_diff.*scaling(ind)./sum(sum(scaling(ind)));
if (~isempty(find(fnew<0)))
    
    neg_ind    = find(fnew < 0);            % finding negative grid points
    new_ind2   = [new_ind, neg_ind'];       % excluding them from mass 
                                            % redistribution as well as
                                            % newly added bins
    ind2       = setdiff([1:N^2],new_ind2);
    mass_rd    = sum(fnew(neg_ind))/length(ind2);
    fnew(ind2) = fnew(ind2) + mass_rd;
    fnew(neg_ind) = 0;                      % setting negative bins to zero
    
end 

return

