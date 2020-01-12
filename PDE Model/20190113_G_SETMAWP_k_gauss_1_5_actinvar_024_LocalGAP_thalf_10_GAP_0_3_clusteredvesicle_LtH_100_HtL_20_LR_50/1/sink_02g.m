% [endo_cell.vSNARE(ind), vSNARE] = ...
%           sink_02b_exclusion (endo_cell.vSNARE(ind),endo_cell.posx(ind),...
%           endo_cell.posy(ind),vSNARE,dx,endo_cell.dx(ind),N,dt,Dconst,1);

function [ cellC, C ] = sink_02g (ind, Xi, Yi, interpD, cellC, posx, posy, C, dx, celldx, dt, Dconst, snk )


% cellC = endo_cell concentrations of protein
% C     = membrane distributionof protein

% number of points of endocytosis
nCell = length(posx);

for i=1:nCell       % for all the points of endocytosis
    
    %% interpolation
    Ci = [];
    if interpD(i) == 2
        ind3 = ind(i,:);
        Ci  = reshape(C(ind3),2,2);
        ptC = interp2(Xi((2*i)-1:2*i,:),Yi((2*i)-1:2*i,:),Ci,posx(i),posy(i));
        
    elseif interpD(i) == 1
        ind3 = ind(i,1:2);
        Ci   = C(ind3);
        ptC  = interp1(Xi((2*i)-1,:),Yi((2*i)-1,:),Ci,posx(i),posy(i));
        
    else
        ind3 = ind(i,1);
        ptC = C(posx(i)/dx,posy(i)/dx);
    end
    
    if isnan(ptC)  % need the zero row or column
        Xi((2*i)-1:2*i,:)
        Yi((2*i)-1:2*i,:)
        Ci
        ind3
    end
    
    %% calculating flux
    % if sink update bin protein concentration
    if snk == 0                       % if not sink:
        % to diffuse in and out use first two lines
        % for reverse sink use second two lines
        
        flux     = ptC - cellC(i);    % +ve if vesicle concentration increasing
        cellC(i) = ptC;               % replace old vesicle concentration
        %         flux     = - dt*4*Dconst*cellC(i)/celldx(i)/celldx(i);   % concentration diffusing out of sink
        %         cellC(i) = cellC(i) + flux;                              % update vesicle concentration
        
    else                              % if sink
        flux     = dt*4*Dconst*ptC/celldx(i)/celldx(i);     % concentration diffusing into sink
        cellC(i) = cellC(i) + flux;        % update vesicle concentration
    end
    
    %% redistributing mass
    if interpD(i) == 2
        
        % flux conservation; re-distribute protein into 4 neighboring bins
        dist(ind3) = sqrt((reshape(Xi((2*i)-1:2*i,:)',1,4)-posx(i)).^2 ...
            +  (reshape(Yi((2*i)-1:2*i,:)',1,4)-posy(i)).^2);
        
    elseif interpD(i) == 1
        dist(ind3) = sqrt((reshape(Xi((2*i)-1,:)',1,2)-posx(i)).^2 ...
            +  (reshape(Yi((2*i)-1,:)',1,2)-posy(i)).^2);
        
    else
        dist(ind3) = 1;
    end
    
    dist(ind3) = 1./dist(ind3);
    fluxm     = flux*(celldx(i)/dx)^2;    % scale mass to a membrane concentration
    C(ind3)    = C(ind3) - fluxm.*dist(ind3)./(sum(dist(ind3)));
    
    % adjust if concentration dips below zero
    ind2      = find(C<0);
    negC      = sum(C(ind2));
    C(ind2)   = 0;
    cellC(i)  = cellC(i) + (dx/celldx(i)).^2 *negC;
end
%%

end



