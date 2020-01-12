      % transform all concentrations so vesicle event is in the centre
      if indx > N/2
        vSNARE   = [vSNARE(  indx-49:end,:) ; vSNARE(  1:indx-50,:)];
        Cdc42T   = [Cdc42T(  indx-49:end,:) ; Cdc42T(  1:indx-50,:)];
        BemGEF42 = [BemGEF42(indx-49:end,:) ; BemGEF42(1:indx-50,:)];
        BemGEF   = [BemGEF(  indx-49:end,:) ; BemGEF(  1:indx-50,:)];
        Cdc42D   = [Cdc42D(  indx-49:end,:) ; Cdc42D(  1:indx-50,:)];
      elseif indx < N/2
        vSNARE   = [vSNARE(  51+indx:end,:) ; vSNARE(  1:50+indx,:)];
        Cdc42T   = [Cdc42T(  51+indx:end,:) ; Cdc42T(  1:50+indx,:)];
        BemGEF42 = [BemGEF42(51+indx:end,:) ; BemGEF42(1:50+indx,:)];
        BemGEF   = [BemGEF(  51+indx:end,:) ; BemGEF(  1:50+indx,:)];
        Cdc42D   = [Cdc42D(  51+indx:end,:) ; Cdc42D(  1:50+indx,:)];
      end
      
      if indy > N/2
        vSNARE   = [vSNARE(  :,indy-49:end) , vSNARE(  :,1:indy-50)];
        Cdc42T   = [Cdc42T(  :,indy-49:end) , Cdc42T(  :,1:indy-50)];
        BemGEF42 = [BemGEF42(:,indy-49:end) , BemGEF42(:,1:indy-50)];
        BemGEF   = [BemGEF(  :,indy-49:end) , BemGEF(  :,1:indy-50)];
        Cdc42D   = [Cdc42D(  :,indy-49:end) , Cdc42D(  :,1:indy-50)];
      elseif indy < N/2
        vSNARE   = [vSNARE(  :,51+indy:end) , vSNARE(  :,1:50+indy)];
        Cdc42T   = [Cdc42T(  :,51+indy:end) , Cdc42T(  :,1:50+indy)];
        BemGEF42 = [BemGEF42(:,51+indy:end) , BemGEF42(:,1:50+indy)];
        BemGEF   = [BemGEF(  :,51+indy:end) , BemGEF(  :,1:50+indy)];
        Cdc42D   = [Cdc42D(  :,51+indy:end) , Cdc42D(  :,1:50+indy)];
      end   