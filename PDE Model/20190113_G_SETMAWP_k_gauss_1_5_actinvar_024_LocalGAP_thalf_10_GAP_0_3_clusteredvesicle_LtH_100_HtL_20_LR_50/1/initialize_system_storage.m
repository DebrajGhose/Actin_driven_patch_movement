
%% Initialize storage for the output.
clear 'c42cen' 'recacen' 'k42_store' 'kra_store' 'endo_cell' 'cable' 'probv_endo';
totcnt = 1;
curr_t = 0;
fpt_flag=0; 
tsi=0;
Fwarp = zeros(N+1); % storage for use in interpolation that needs wrap around func values
c42cen=[];recacen=[]; 
endo_cell        = struct('endo_time',[],'pick_time',[],'vSNARE',[],'Cdc42T',[], ...
    'Cdc42D',[], 'BemGEF42',[], 'BemGEF',[], 'posx',[],'posy',[],'dx',[],'num',0);
mass_store    = fopen(['data_mass_run' num2str(runi) '.txt'],                  'w'); % going to hold total mass
Internalcompartment = 0.7;                  % relative size (area) of internal compartment vs plasma membrane 
dx  = Cellsize*sqrt(pi)/N;                  % spatial meshsize (1D)
dx2 = sqrt(Internalcompartment)*dx;
eta0=0.01;                     %eta: Vm/Vc, membrane/cytoplasm volume correction
eta =eta0;                     % notes in '08_endocytic_vesicle_30Dec2010.doc'
mem_depth  = (Cellsize/2)*(1-(1/((eta0+1)^(1/3))));
cyt_mult   = ((Cellsize/2)^3)/(eta0+1);
Rnew_mult  = N/2/(pi^0.5);

%% initial concentrations

load('Init_AllProteins_Pher_10nM.mat');

%% Calclate the global mass __________________________________________
Cdc42s = Cdc42Dc ...                                    % Cdc42D cytoplasmic content 
+ eta*(mean(mean(Cdc42D+Cdc42T+BemGEF42))) ...   % Cdc42 membrane content 
+ eta*(Cdc42Dic*(dx2^2)/(dx^2)) ...                    % Cdc42D on inner membrane           
+ eta*sum(endo_cell.Cdc42T(1:endo_cell.num)   .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2) ...     % Cdc42T in vesicle                            
+ eta*sum(endo_cell.Cdc42D(1:endo_cell.num)   .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2) ...     % Cdc42D in vesicle 
+ eta*sum(endo_cell.BemGEF42(1:endo_cell.num) .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2);        % BemGEF42 in vesicle 

Bem1s = BemGEFc ...                                    % cytoplasmic content scaled to outer membrane
+ eta*(mean(mean(BemGEF+BemGEF42))) ...                % Bem1 membrane content  
+ eta*sum(endo_cell.BemGEF(1:endo_cell.num)   .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2) ...     % Bem1 vesicle content 
+ eta*sum(endo_cell.BemGEF42(1:endo_cell.num) .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2);        % BemGEF42 vesicle content                                                        

vSNAREs        = eta*mean(mean(vSNARE)) ...              % v-SNARE membrane content
+ eta*(vSNAREic*(dx2^2)/(dx^2)) ...                   % v-SNARE on inner membrane content 
+ eta*sum(endo_cell.vSNARE(1:endo_cell.num) .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2);          % v-SNARE vesicle content 

%% Clean vesicle containers endo and exo
cable       = [];       % hold upto MAXcable cable co-ordinates
% record_vSNARE = [];     % record [v-SNARE] at time of endocytosis
% record_endot = [];      % record time of endocytosis
% record_lifetime = [];   % record lifetime of each sink
% record_Rec = [];     % record [Rec] at time of endocytosis
% record_RecA = [];     % record pheromone bound receptor [RecA] at time of endocytosis

%exocytic thresholds for hill distribution
if(max(max(RecA))~=0)     k1=0.5*max(max(RecA)); else    k1=1; end %exponent and threshold for RecA, threshold shosen as a half of <Rec> in the absense of pheromone 1.5 microMolar
k2=0.5*max(max(Cdc42T+BemGEF42)); %about a 1/2 of the maximum Cdc42T concentration in a well formed peak in our simulations
k42_store=zeros([1 ceil(sim/tsave)+1]);
k42_store(1)=k2;
% kra_store=zeros([1 ceil(sim/tsave)+1]);
% kra_store(1)=k1;
if(max(max(vSNARE))~=0)     kv=0.5*max(max(vSNARE)); else    kv=1; end 
kvsn_store=zeros([1 ceil(sim/tsave)+1]);
kvsn_store(1)=kv;

%% report at the end
display(['--> Initialization lasted ' num2str(toc(tot_runi_start)) ' (s). Starting simulations...']); %Masha debug
