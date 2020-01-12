Adj=[];
corrcoef=1.3;
%% adjustible step size with Non-negative concentrations control
% Reminder: nRsteps = 100;      % take nRsteps reaction steps per diffusion step
% REMINDER: dt2 = dt/nRsteps;   % reaction time-step
% tStart3  = tic;
react_t=0;
adj_count=0;
while (react_t<=dt)
    dts_i=dt2;
    
    Reaction_Euler_step;
    
   % Species_h  = Species    + dts_i*react_flux(dts,Species);
    while (~isempty(find(Cdc42Th<0))  | ~isempty(find(BemGEF42h<0)) | ~isempty(find(BemGEFh<0)) ...
          | ~isempty(find(BemGEFch<0)) | ~isempty(find(Cdc42Dh<0))   | ...
             ~isempty(find(Cdc42Dch<0))  ...
          | ~isempty(find(Cdc42Dich<0)) )
        dts_i=dts_i/1.5;
        adj_count=adj_count+1; %Masha debug
    
        Reaction_Euler_step;
        
    end     
    react_t=react_t+dts_i;
    Cdc42T   = Cdc42Th;
    BemGEF42 = BemGEF42h;
    BemGEF   = BemGEFh;
    BemGEFc  = BemGEFch;
    Cdc42D   = Cdc42Dh;
    Cdc42Dc   = Cdc42Dch;
    Cdc42Dic = Cdc42Dich;
   
end 
Adj=[Adj adj_count];