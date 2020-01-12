

%Take current profile to estimate HIll Coefficients

k2=0.5*max(max(Cdc42T+BemGEF42)); %about a 1/2 of the maximum Cdc42T concentration in a well formed peak in our simulations
if(max(max(vSNARE))~=0),     kv=0.5*max(max(vSNARE)); else    kv=1; end


% cargo "sinks" into cells destined for endocytosis

if ( endo_cell.num > 0 )
    %  endo_time includes that 15-s endo_red_time
    %  find "green" bins
    ind = find(endo_cell.endo_time < 0);
    
    [ind2,Xi,Yi,interpD] = sink_ind_02g (endo_cell.posx(ind),endo_cell.posy(ind), dx, N);
    
    [endo_cell.vSNARE(ind), vSNARE] = ...
        sink_02g (ind2,Xi,Yi,interpD, ...
        endo_cell.vSNARE(ind),endo_cell.posx(ind),...
        endo_cell.posy(ind),vSNARE,dx,endo_cell.dx(ind),dt,Dconst,1);
    
    [endo_cell.Cdc42T(ind), Cdc42T] = ...
        sink_02g (ind2,Xi,Yi,interpD, ...
        endo_cell.Cdc42T(ind),endo_cell.posx(ind),...
        endo_cell.posy(ind),Cdc42T,dx,endo_cell.dx(ind),dt,Dconst,0);
    [endo_cell.Cdc42D(ind), Cdc42D] = ...
        sink_02g (ind2,Xi,Yi,interpD, ...
        endo_cell.Cdc42D(ind),endo_cell.posx(ind),...
        endo_cell.posy(ind),Cdc42D,dx,endo_cell.dx(ind),dt,Dconst,0);
    [endo_cell.BemGEF42(ind),  BemGEF42] = ...
        sink_02g (ind2,Xi,Yi,interpD, ...
        endo_cell.BemGEF42(ind),endo_cell.posx(ind),...
        endo_cell.posy(ind),BemGEF42,dx,endo_cell.dx(ind),dt,Dconst,0);
    
    [endo_cell.BemGEF(ind), BemGEF] = ...
        sink_02g (ind2,Xi,Yi,interpD, ...
        endo_cell.BemGEF(ind),endo_cell.posx(ind),...
        endo_cell.posy(ind),BemGEF,dx,endo_cell.dx(ind),dt,Dconst,0);
    
end
ticd=tic;
% diffuse
% membrane
vSNARE   = reshape(Hopvsnare\reshape(vSNARE,  N^2,1),N,N);
Cdc42T   = reshape(Hop\reshape(Cdc42T,  N^2,1),N,N);
BemGEF42 = reshape(Hop\reshape(BemGEF42,N^2,1),N,N);
BemGEF   = reshape(Hop\reshape(BemGEF,  N^2,1),N,N);
Cdc42D   = reshape(Hop\reshape(Cdc42D,  N^2,1),N,N);

k2bvesc   = reshape(Hop\reshape(k2bvesc,  N^2,1),N,N); %this is being stretched and pinched but not geting endocytosed

% calculate new cytoplasmic concns by averaging over all cells
BemGEFc(:,:) = mean(mean(BemGEFc));
Cdc42Dc(:,:)  = mean(mean(Cdc42Dc));

%tdiffuse=[tdiffuse toc(ticd)];
%% exocytosis
%Find cable attachment zones, calculate probabilities
% Cable attach
if length(cable)/2 < MAXcable                           % attach new cables
    cable_add = [];                                    % store new cable positions
    for  i = 1:(MAXcable - (length(cable)/2))
        if random('unif',0,1) <= attachPROB           % attach a cable
            %smooth overlap for cable attachments + exclude the locations with endocytic patches
            ct=Cdc42T+BemGEF42;
            p42=(ct.^n2)./(k2^n2 + ct.^n2); %probability to find Cdc42T modeled via Hill function
            pcab=p42;
            
            %Actin may increase probability of binding cables near them.
            %Introduce Gaussian around cable centroid.
            
            if numel(cable)>0 %cables HAVE to be present for Gaussian to have any effect
               
                [~,max42ind]=max(Cdc42T(:)); [max42indr,max42indc] = ind2sub([N,N],max42ind); %find position of max 42 because that's easier/faster than finding centroid for every dt
                
                cablecen = calculate_cable_centroid(cable,[max42indr,max42indc],N); %find centroid of cables
                modactrow = mod(actrowcoord + 50 - cablecen(end,1) ,  N+10^-12 );
                modactcol = mod(actcolcoord + 50 - cablecen(end,2), N+10^-12 );
                Pact = k_gauss*exp( -( (modactrow-50).^2/(2*actinsd^2) + (modactcol-50).^2/(2*actinsd^2) ) ); %generate a gaussain where you think actin cable centroid is
                
            else, Pact = 0; %if no cables, they should have no effect. This needs to be a matrix to avoid dimension mismatch error.
            end
            
            pcab = pcab + pcab.*Pact; %actin cable centroid affects probability of attachment
            
            %remove any possibility to attach a cable at the existing endocytic patch
            if ( endo_cell.num > 0 )
                for i=1:endo_cell.num
                    pcab(ceil(endo_cell.posx(i)./dx),ceil(endo_cell.posy(i)./dx))=0;
                end
            end
            %remove any possibility of attaching 2 cables on top of each other
            indcab=1:2:(length(cable)-1);
            pcab(sub2ind(size(pcab),cable(indcab),  cable(indcab+1)   ))=0;
            pcab=pcab/sum(sum(pcab));
            pcabl=zeros([1 N^2]);sumi=0;
            %inds=randperm(N^2); %shuffle matrix indexes to calc cumul. prob in random order
            inds=1:1:N^2;
            for (i=inds)
                sumi=sumi+pcab(i);
                pcabl(i)=sumi;
            end
            dice3=random('unif',0,1,1,1);
            for (k=1:length(pcabl))
                if (pcabl(k)>dice3)
                    subcab=k;
                    break;
                end
            end
            [indx,indy]=ind2sub(N,inds(subcab));
            
            cable_add = [cable_add, indx, indy];      % store new cable positions
        end
    end
end

% Cable detatch
if length(cable) ~= 0                                  % detatch cables
    cable_rmv = [];
    for  i = 1:2:length(cable)
        if random('unif',0,1) <= detatchPROB            % detatch a cable
            cable_rmv = [cable_rmv,i,i+1];
        end
    end
    cable(cable_rmv) = [];
end

cable = [cable cable_add];
cable_add = [];

%% an exocytotic event may happen

dice = random('unif',0,1); % dice=random number between 0 and 1
if ( dice < dt/periodIn && curr_t && length(cable)>0) % exocytosis will happen if dice<dt/periodIn (=0.05/2.4=0.02: 2% chance)
    exo_cnt = exo_cnt + 1;
    
    % Pick cable
    dice = random('unif',0,1,1,(length(cable)));       %Used to be length(cable)/2, but corrected it -- Debraj
    i    = find(dice==max(dice));
    if (mod(i,2)==0) %if the index is 2th in a found pair (Y coordinate of a cable), shift to the X index.
        i=i-1;
    end
    
    % Randomize the cable location for vesicle fusion within 4  4pnt combinations around the original location
    quadi=randi(4,1); % randomly choose a quadrant notrthwest-1, northeast -2, southwest -3, southeast -4.
    switch quadi %depending on a chosen quadrant always put a (zero)-fusion point to be the southwest corner of the 4 grid points. fusion will happen to the upper east
        case 1
            indx= mod(cable(i)-1,N);
            if (indx==0) indx=N; end
            indy=mod(cable(i+1),N);
            if (indy==0) indy=N; end
        case 2
            indx= mod(cable(i),N);
            if (indx==0) indx=N; end
            indy=mod(cable(i+1),N);
            if (indy==0) indy=N; end
        case 3
            indx= mod(cable(i),N);
            if (indx==0) indx=N; end
            indy=mod(cable(i+1)-1,N);
            if (indy==0) indy=N; end
        case 4
            indx= mod(cable(i)-1,N);
            if (indx==0) indx=N; end
            indy=mod(cable(i+1)-1,N);
            if (indy==0) indy=N; end
        otherwise
            display(['quadi = ' num2str(quadi) ': could not choose the right quandrant to fuse a vesicle']);
    end
    indxi=mod(([indx indx indx indx] + [0 0 1 1]),N); %indexes for a north-east 3 points from the original fusion location
    indyi=mod(([indy indy indy indy] + [0 1 1 0]),N); %indexes for a north-east 3 points from the original fusion location. All interpolations are hard wired to fit, stretch and approximate those 3upper-right pnts
    indxi(find(indxi==0))=N;
    indyi(find(indyi==0))=N;
    inds=sub2ind(size(vSNARE), indxi, indyi); %subscripts of grid points in which to insert vesicle content
    % Radially 'stretch' [Cdc42] distribution, set up variables for interpolation
    
    r_exo  = nVInlen*dx/(pi^0.5);  % radius of area added - old membrane size
    % area added is (nVInlen*dx)^2
    C = N*dx;
    hole_radial_exo_area_05;                      % call script 'hole' to compute
    
    % transform all concentrations so vesicle event is in the centre
    vesicle_position_center;
    
    % add naked membrane to external membrane
    vSNARE   = exo_stretch_radial_01 ( vSNARE,  Fwarp,  x0, y0, x1, y1,[N/2,(N/2)+1],[N/2,(N/2)+1],N^2/(N^2+nVInlen^2));
    Cdc42T   = exo_stretch_radial_01 ( Cdc42T,  Fwarp,  x0, y0, x1, y1,[N/2,(N/2)+1],[N/2,(N/2)+1],N^2/(N^2+nVInlen^2));
    BemGEF42 = exo_stretch_radial_01 ( BemGEF42,Fwarp,  x0, y0, x1, y1,[N/2,(N/2)+1],[N/2,(N/2)+1],N^2/(N^2+nVInlen^2));
    BemGEF   = exo_stretch_radial_01 ( BemGEF,  Fwarp,  x0, y0, x1, y1,[N/2,(N/2)+1],[N/2,(N/2)+1],N^2/(N^2+nVInlen^2));
    Cdc42D   = exo_stretch_radial_01 ( Cdc42D,  Fwarp,  x0, y0, x1, y1,[N/2,(N/2)+1],[N/2,(N/2)+1],N^2/(N^2+nVInlen^2));
    
    k2bvesc   = exo_stretch_radial_01 ( k2bvesc,  Fwarp,  x0, y0, x1, y1,[N/2,(N/2)+1],[N/2,(N/2)+1],N^2/(N^2+nVInlen^2)); % %this is being stretched and pinched but not geting endocytosed
    
    clear xold yold
    
    % transform all concentrations back
    vesicle_position_restore;
    
    %% move endocytic patches after exocytosis
    dx_new  = dx*sqrt(N^2+nVInlen^2)/N;
    C = N*dx_new;                  % move patches on new membrane
    endo_cell = moving_patches_exo_radial(r_exo,(indx*dx_new)+(dx_new/2),(indy*dx_new)+(dx_new/2),endo_cell,C);
    
    %% adjust size of external membrane, and diffusion operator
    dx  = dx*sqrt(N^2+nVInlen^2)/N;
    eta = ((Rnew_mult^3 * dx^3) - (((Rnew_mult * dx)-mem_depth)^3))/cyt_mult;
    Dxx = Dxx*N^2/(N^2+nVInlen^2);
    Hop = speye(N^2) - kron(speye(N),Dxx) - kron(Dxx,speye(N));
    
    
    %% do the same with vsnare diffusion operator
    
    Dxxvsnare = Dxxvsnare*N^2/(N^2+nVInlen^2);
    Hopvsnare = speye(N^2) - kron(speye(N),Dxxvsnare) - kron(Dxxvsnare,speye(N));
    
    
    %% adjust size of internal membrane
    dx2new = dx2*sqrt(N^2-(nVInlen*dx/dx2)^2)/N;
    if (indx<=0 | indx > N | indy <=0 | indy>N)
        display(['dice= ' num2str(dice) 'indx = ' num2str(indx) ' indy = ' num2str(indy)]);
    end
    
    %% Set [cargo] in new patch of membrane
    %  NEw june 2014 new way
    vSNARE(inds) = 10*vSNAREic;  %pbc wrap
    Cdc42D(inds) = Cdc42Dic;
    k2bvesc(inds) = k2bdeliv; % Add GAP activity to delivery location
    
    %% remove transferred cargo from internal compartment
    vSNAREic = vSNAREic*(N^2*dx2^2-40*dx^2)/(N^2*dx2new^2);
    dx2 = dx2new;
    
end

%% an endocytotic event may happen
%Edited to assure proper spatial distribution of endo patches
dice = random('unif',0,1); % dice=random number between 0 and 1
if ( dice(1) < dt/periodOut && curr_t ) % endocytosis will happen if dice<dt/periodOut (=0.05/0.6=0.08: 8% chance)
    probv_endo=vSNARE.^nv./(kv^nv+vSNARE.^nv);
    probv_endo(find(probv_endo<max(max(probv_endo))/40))=max(max(probv_endo))/40; % maximum ratio is set to 40, to still have some posibility of endocytosis in the back
    %remove any possibility to attach a cable at the existing endocytic patch
    if ( endo_cell.num > 0 )
        %ind = find(endo_cell.endo_time < 0);
        for i=1:endo_cell.num
            probv_endo(ceil(endo_cell.posx(i)./dx),ceil(endo_cell.posy(i)./dx))=0;
        end
    end
    %remove any possibility of attaching 2 cables on top of each other
    if(length(cable)>1)
        indcab=1:2:(length(cable)-1);
        probv_endo(sub2ind(size(probv_endo),cable(indcab),         cable(indcab+1)          ))=0;
        probv_endo(sub2ind(size(probv_endo),mod(cable(indcab),N)+1,cable(indcab+1)          ))=0;
        probv_endo(sub2ind(size(probv_endo),cable(indcab),         mod(cable(indcab+1),N)+1 ))=0;
        probv_endo(sub2ind(size(probv_endo),mod(cable(indcab),N)+1,mod(cable(indcab+1),N)+1 ))=0;
    end
    probv_endo=probv_endo/sum(sum(probv_endo));
    pendol=zeros([1 N^2]);sumi=0;
    for (i=1:N^2)
        sumi=sumi+probv_endo(i);
        pendol(i) = sumi;
    end
    dice3=random('unif',0,1,1,1);
    for (k=1:length(pendol))
        if (pendol(k)>dice3)
            subcab=k;
            break;
        end
    end
    [indx,indy]=ind2sub(N,subcab);
    
    
    % put this bin onto array with bins destined to endocytose
    % set negative endocytosis time for bins not yet full
    endo_cell.endo_time = [endo_cell.endo_time,-1];
    endo_cell.pick_time = [endo_cell.pick_time,curr_t];
    
    % adding proteins to vesicle
    endo_cell.vSNARE  = [endo_cell.vSNARE,   vSNARE(indx,indy)];
    endo_cell.Cdc42T  = [endo_cell.Cdc42T,   Cdc42T(indx,indy)];
    endo_cell.Cdc42D  = [endo_cell.Cdc42D,   Cdc42D(indx,indy)];
    endo_cell.BemGEF42= [endo_cell.BemGEF42, BemGEF42(indx,indy)];
    endo_cell.BemGEF  = [endo_cell.BemGEF,   BemGEF(indx,indy)];
    
    % recording the position of vesicle
    endo_cell.posx = [endo_cell.posx,indx*dx];
    endo_cell.posy = [endo_cell.posy,indy*dx];
    % recording the dx value at time vesicle was created
    endo_cell.dx   = [endo_cell.dx,dx];
    endo_cell.num  = endo_cell.num+1;
    
    
    % Radially 'squeeze' [Cdc42] distribution, set up variables for interpolation
    %
    
    r_endo  = nVOutlen*dx/(pi^0.5);  % radius of area removed - old membrane size
    % area removed is (nVOutlen*dx)^2
    C = N*dx;
    hole_radial_endo_area_04;                         % call script 'hole' to compute new grid
    % in this case the new grid has a point
    % missing at the endocytosis point
    
    % transform all concentrations so vesicle event is in the centre
    vesicle_position_center;
    
    vSNARE   = endo_squeeze_radial ( vSNARE,   Fwarp,  x0, y0, x1, y1, N/2, N/2, N^2/(N^2-nVOutlen^2));
    Cdc42T   = endo_squeeze_radial ( Cdc42T,   Fwarp,  x0, y0, x1, y1, N/2, N/2, N^2/(N^2-nVOutlen^2));
    BemGEF42 = endo_squeeze_radial ( BemGEF42, Fwarp,  x0, y0, x1, y1, N/2, N/2, N^2/(N^2-nVOutlen^2));
    BemGEF   = endo_squeeze_radial ( BemGEF,   Fwarp,  x0, y0, x1, y1, N/2, N/2, N^2/(N^2-nVOutlen^2));
    Cdc42D   = endo_squeeze_radial ( Cdc42D,   Fwarp,  x0, y0, x1, y1, N/2, N/2, N^2/(N^2-nVOutlen^2));
    
    k2bvesc   = endo_squeeze_radial ( k2bvesc,   Fwarp,  x0, y0, x1, y1, N/2, N/2, N^2/(N^2-nVOutlen^2)); %this is being stretched and pinched but not geting endocytosed
    
    clear xold yold
    
    % transform all concentrations back
    vesicle_position_restore;
    
    %% move endocytic patches after endocytosis
    
    dx_new  = dx*sqrt(N^2-nVOutlen^2)/N;
    C = N*dx_new;                    % move patches on new membrane size
    endo_cell = moving_patches_endo_radial(r_endo,indx*dx_new,indy*dx_new,endo_cell,C);
    
    %%
    
    %% adjust size of membrane, and diffusion operator
    dx  = dx*sqrt(N^2-nVOutlen^2)/N;
    eta = ((Rnew_mult^3 * dx^3) - (((Rnew_mult * dx)-mem_depth)^3))/cyt_mult;
    Dxx = Dxx*N^2/(N^2-nVOutlen^2);
    Hop = speye(N^2) - kron(speye(N),Dxx) - kron(Dxx,speye(N));
    
   
    %% do the same for vnsare diffusion
    Dxxvsnare = Dxxvsnare*N^2/(N^2-nVOutlen^2);
    Hopvsnare = speye(N^2) - kron(speye(N),Dxxvsnare) - kron(Dxxvsnare,speye(N));
    
end

% look for bins that are full, then have them wait 15 s
ind = find( ...
    (endo_cell.vSNARE > patch_fill | curr_t-endo_cell.pick_time > patch_fill_give_up) ...
    & endo_cell.endo_time <= 0);
endo_cell.endo_time(ind) = curr_t + endo_red_time;

% look for bins that need to endocytose now
ind = find((endo_cell.vSNARE>patch_fill | curr_t-endo_cell.pick_time ...
    > patch_fill_give_up) & endo_cell.endo_time<=curr_t+0.5*dt & endo_cell.endo_time>=0);
% if last condition false, that means bins need to turn red and wait 15 s

if (~isempty(ind))  % a cell need to endocytose!
    % length(ind)
    %endotimes=[endotimes curr_t];
    
    % only adding one vesicle to inner membrane despite the face that
    % length(ind) could be greater than one
    mydx2 = sqrt( dx2^2 + sum((nVOutlen*endo_cell.dx(ind)).^2 ) /N^2 );
    
    %      fprintf(rec_count_time,'%5.4f %5.4f %5.4f %5.4f %5.4f\n',mean(mean(Rec)), sum(endo_cell.Rec(ind)), Recic, mean(mean(RecA)),sum(endo_cell.RecA(ind)) );
    
    
    % dump sink content to internal compartment
    vSNAREic = ((vSNAREic*N^2*dx2^2) + ...
        (sum(endo_cell.vSNARE(ind).* (endo_cell.dx(ind).^2) )) )/mydx2^2/N^2;
    
    Cdc42Dic = ((Cdc42Dic*N^2)*dx2^2 + ...
        sum( endo_cell.Cdc42D(ind)                            .*((endo_cell.dx(ind).^2)) ) ...
        + sum( (endo_cell.Cdc42T(ind) + endo_cell.BemGEF42(ind)).*((endo_cell.dx(ind).^2)) ) ...
        )/mydx2^2/N^2;
   
    % adding total Bem1 to cytoplasm
    BemGEFc = BemGEFc + ...
        (eta*sum( (endo_cell.BemGEF42(ind) + endo_cell.BemGEF(ind)).*(endo_cell.dx(ind).^2)/(dx^2) )/(N^2));
    
    dx2 = mydx2;
    
    % take this bin off the array
    
    newind = setdiff([1:endo_cell.num],ind); % values in '[1:endo_cell.num]' = 1,2,3,...,endo_cell.num
    % and not in 'ind'
    endo_cell.endo_time = endo_cell.endo_time(newind);
    endo_cell.pick_time = endo_cell.pick_time(newind);
    
    endo_cell.vSNARE    = endo_cell.vSNARE(newind);
    endo_cell.Cdc42T    = endo_cell.Cdc42T(newind);
    endo_cell.Cdc42D    = endo_cell.Cdc42D(newind);
    endo_cell.BemGEF42  = endo_cell.BemGEF42(newind);
    endo_cell.BemGEF    = endo_cell.BemGEF(newind);
    
    endo_cell.posx      = endo_cell.posx(newind);
    endo_cell.posy      = endo_cell.posy(newind);
    endo_cell.dx        = endo_cell.dx(newind);
    endo_cell.num       = endo_cell.num - length(ind);
    
end

