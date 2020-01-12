      % transform all concentrations back
      if indx > N/2
        vSNARE   = [vSNARE(  151-indx:end,:) ; vSNARE(  1:150-indx,:)];
        Cdc42T   = [Cdc42T(  151-indx:end,:) ; Cdc42T(  1:150-indx,:)];
        BemGEF42 = [BemGEF42(151-indx:end,:) ; BemGEF42(1:150-indx,:)];
        BemGEF   = [BemGEF(  151-indx:end,:) ; BemGEF(  1:150-indx,:)];
        Cdc42D   = [Cdc42D(  151-indx:end,:) ; Cdc42D(  1:150-indx,:)];
      elseif indx < N/2
        vSNARE   = [vSNARE(  51-indx:end,:) ; vSNARE(  1:50-indx,:)];
        Cdc42T   = [Cdc42T(  51-indx:end,:) ; Cdc42T(  1:50-indx,:)];
        BemGEF42 = [BemGEF42(51-indx:end,:) ; BemGEF42(1:50-indx,:)];
        BemGEF   = [BemGEF(  51-indx:end,:) ; BemGEF(  1:50-indx,:)];
        Cdc42D   = [Cdc42D(  51-indx:end,:) ; Cdc42D(  1:50-indx,:)];
      end
      
      if indy > N/2
        vSNARE   = [vSNARE(  :,151-indy:end) , vSNARE(  :,1:150-indy)];
        Cdc42T   = [Cdc42T(  :,151-indy:end) , Cdc42T(  :,1:150-indy)];
        BemGEF42 = [BemGEF42(:,151-indy:end) , BemGEF42(:,1:150-indy)];
        BemGEF   = [BemGEF(  :,151-indy:end) , BemGEF(  :,1:150-indy)];
        Cdc42D   = [Cdc42D(  :,151-indy:end) , Cdc42D(  :,1:150-indy)];
      elseif indy < N/2
        vSNARE   = [vSNARE(  :,(51-indy):end) , vSNARE(  :,1:(50-indy))];
        Cdc42T   = [Cdc42T(  :,(51-indy):end) , Cdc42T(  :,1:(50-indy))];
        BemGEF42 = [BemGEF42(:,(51-indy):end) , BemGEF42(:,1:(50-indy))];
        BemGEF   = [BemGEF(  :,(51-indy):end) , BemGEF(  :,1:(50-indy))];
        Cdc42D   = [Cdc42D(  :,(51-indy):end) , Cdc42D(  :,1:(50-indy))];
      end    