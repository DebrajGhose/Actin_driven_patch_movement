% Euler_step

% membrane
Cdc42Th    = Cdc42T    + dts_i.* ((k2a*BemGEF+k3*BemGEF42).*Cdc42D - k2btotal.*Cdc42T - k4a*BemGEF.*Cdc42T + k4b*BemGEF42 - k7*BemGEFc.*Cdc42T);

BemGEF42h     = BemGEF42     + dts_i.*(k4a*BemGEF.*Cdc42T - k4b*BemGEF42 + k7*BemGEFc.*Cdc42T);             

BemGEFh    = BemGEF    + dts_i.*(k1a*BemGEFc - k1b*BemGEF - k4a*BemGEF.*Cdc42T + k4b*BemGEF42);

Cdc42Dh    = Cdc42D    + dts_i.*(k5a*Cdc42Dc - k5b*Cdc42D + k2btotal.*Cdc42T - (k2a*BemGEF+k3*BemGEF42).*Cdc42D );

% cytoplasm

BemGEFch    = BemGEFc    + dts_i.*mean(mean((eta*(k1b*BemGEF - (k1a+k7*Cdc42T).*BemGEFc))));

Cdc42Dch  = Cdc42Dc  + dts_i.*mean(mean((eta*(k5b*Cdc42D  - k5a*Cdc42Dc) + ...
                     eta*(k5b*Cdc42Dic - k5a*Cdc42Dc   )*dx2*dx2/dx/dx )));

%internal compartment

Cdc42Dich = Cdc42Dic + dts_i.*(k5a*mean(mean(Cdc42Dc)) - k5b*Cdc42Dic);
