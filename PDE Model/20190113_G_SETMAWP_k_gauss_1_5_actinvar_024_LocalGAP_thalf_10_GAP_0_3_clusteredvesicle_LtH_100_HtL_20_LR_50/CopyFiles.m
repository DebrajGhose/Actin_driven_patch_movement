%create multiple files for simulation

rng(753655) %seed for random generator 

sims = 400; %number of simulations
tic
for i = 2:sims


copyfile('1',num2str(i)) ;


end

%generate random seeds for all the files

    for j = 1:sims
    
    seed = randi([0,2^32-1]);
    
    %save the seed file
    
    
    filename = strcat(num2str(j),'/seed.mat');
   
    save(filename,'seed');
    
    
    end
toc