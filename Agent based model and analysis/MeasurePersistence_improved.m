%% Analyze patch step characteristics and persistence
clear all
close all
sprintf(['\n']) %create some white space for better readability
clean_data = 0; %Go through cleaning data process
generate_in_silico_data = 0; %Collect persistence value from insilico model and simulate on a sphere
plot_persistence = 0; %Plot persistence measurements
persistenceautocorrelation = 0; %see what autocoreelation graph looks like for persistent movement
persistenceaudocorrelation_experimental_data = 0;
plot_persistence_2 = 1; % Plot step length and persistence, but with preloaded data
plot_persistence_3 = 0; % Plot step length, persistence, and a summary of your results
difference_between_polarisome_and_patch = 0; %find distance between actin cable centroid and your patch
plot_diffusion_constants = 0; %plot diffusion constants for the experiments and simulations
simulate_search_efficiency = 0; %see how efficient different search methods are
simulate_search_efficiency_gamma_dist = 0; %see how efficient different search methods are while sampling from different gamma distributions
create_figure_of_cells_searching = 0; %create a figure for particles moving along surface of a sphere
see_patch_variability = 0; %see how patch height/amount of Cdc42 changes for different parameters
model_microscope_point_spread = 0; %model how point spread function might affect visibility of actin cables
sample_profiles_of_polarisome_and_patch = 0;
make_figure_bursty_deliv = 0;

%% Set up colors


newcolors = [ 203 91 90;...
    0 0 0;...
    83 187 110;...
    130 172 162;...
    76 180 209;...
    45 112 178;
    26 110 38]/255; 


%{

newcolors = [ 225 151 76 ;... %orange
    0 0 0 ;... %black
    114 147 203 ;... %blue
    132 186 91 ;... %green
    144 103 167 ; ... %purple
    171 104 87 ;... %maroon
    211 94 96 ;... %red
    128 133 133;...%gray
    204 194 16]/255; %gold

%}

set(groot,'defaultAxesColorOrder',newcolors)

%% Clean data

Datafrom = 'Data from Allie 150921_latA';
pherandgenotype = 'lata';

if clean_data == 1
    
    alltracks = {}; %store your tracks here
    origintrack = readtable( [Datafrom,'\',pherandgenotype,'\origin.csv'] );
    origintrack = [ cellfun(@str2num,origintrack{1:end,5}) , cellfun(@str2num,origintrack{1:end,18:20}) ] ; %extract timepoints and the xyz positions
    
    for cellno = 1:64
        
        if exist([Datafrom,'\',pherandgenotype,'\cell',num2str(cellno),'.csv'],'file')
            celltrack = readtable( [Datafrom,'\',pherandgenotype,'\cell',num2str(cellno),'.csv'] );
            celltrack = [ cellfun(@str2num,celltrack{1:end,5}) , cellfun(@str2num,celltrack{1:end,18:20}) ] ; %extract timepoints and the xyz positions
            diffa = diff(celltrack(:,1)); %throw out timepoints where there are two patches
            diffa(end+1) = diffa(end);
            elimtheserep = diffa==0; %throw awaytimepoints where there are multiple patches
            elimtheserep2 = logical([elimtheserep(1) ; elimtheserep(1:end-1)]);
            elimtheserep = elimtheserep|elimtheserep2;
            celltrack(elimtheserep,:) = [];
            
            elimtimepoints = setdiff( origintrack(:,1),celltrack(:,1) ); %find timepoints that are missing
            originalt = origintrack; %make copy of origin matrix
            originalt(elimtimepoints,:) = []; %remove missing timepoints from origin matrix
            celltrack(:,2:4) = celltrack(:,2:4)-originalt(:,2:4); %subtract coordinates by origin to stabilize patch locations
            
            %address discontinuities in path
            
            while size(celltrack,1)>0 %I am going through the matrix and looking for discontinuities in the timepoint order.
                %each time I find a discontinuity, I move elements upstream of the
                %discontinuity into a separate cell matrix. I stop performing this
                %operation once there are no elements left in a.
                
                diffa = diff(celltrack(:,1))'; %if patch disappears at some timepoint, you want to end the track right there
                
                indices = [find(diffa>1) , numel(celltrack(:,1))]'; %find indices you want to eliminate rows up to.
                %In case there are no discontinuities, indices will simply contain,
                %the row lenght of the matrix and copy the entire thing into a cell
                %matrix.
                
                if size(celltrack(1:indices(1),:),1)~=1 %don't bother storing tracks that are only one element long
                    
                    alltracks{end+1} = celltrack(1:indices(1),:); %add another track to the cell matrix
                end
                celltrack(1:indices(1),:) = []; %remove track that has already been stored
                
            end
        end
    end
    
    save(['AllTracks_', Datafrom , pherandgenotype],'alltracks');
    
end

%% Generate pesistence data from in silico model and plot

if generate_in_silico_data == 1
    
    foldername ='C:\Users\dg144\Desktop\Projects\Simulations\20190113_SETMAWP_Treadmill\OUT\' ;
    for filenamecell = {'20190113_G_SETMAWP_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_clusteredvesicle_LtH_100_HtL_20_LR_100_k3_06645'}
        
        filename = char(filenamecell);
        if ~strcmp(filename,'SampleFromUniformAngle'), load([foldername,filename,'.mat']); end
        
        skip = 18; %make time interval 90s
        cellradius = 2.5; %microns
        insilicopixtomic = 0.0886;
        
        if strcmp(filename,'SampleFromUniformAngle')
            
            load('InVivo.mat','allrs');
            vectr = allrs;
            
        else
            
            allx = allx(:,1:skip:end)*insilicopixtomic; ally = ally(:,1:skip:end)*insilicopixtomic;
            vectx = allx(:,2:end) - allx(:,1:end-1); vecty = ally(:,2:end) - ally(:,1:end-1); % find vectors
            [vecttheta,vectr] = cart2pol(vectx,vecty);
            diffangs = (angdiff(vecttheta(:,2:end),vecttheta(:,1:end-1)));
            
        end
        
        %diffangs = normrnd(0,0.001,[400,22]); %for testing and debugging
        %now use these angles and rs to simulate a diffusing particle on a sphere
        
        numbsteps = 30000; %number of simulations steps you want this to take
        storedpath = zeros(numbsteps,3);
        
        storedpath(1,:) = [ 0 cellradius 0 ];
        
        for iindex = 2:numbsteps
            
            if strcmp(filename,'SampleFromUniformAngle')
                ang = rand()*pi;
            else
                ang = randsample(diffangs(:),1);
            end
            
            mydist = randsample(vectr(:),1);
            
            [storedpath] = Function_sphere_walk( mydist , ang , storedpath , iindex , cellradius);
            
            if 0 %choose whether to draw traces
                %draw things
                hold off
                
                if ~exist('videowrite','var')
                    videowrite = VideoWriter('persistentmotionexample.avi');
                    open(videowrite);
                end
                
                [xds,yds,zds]  = sphere;
                surf(  0.99*cellradius*xds ,  0.99*cellradius*yds ,  0.99*cellradius*zds  ,'FaceColor' , [ 0.9 0.7 0.9 ] , 'EdgeColor' , [1 1 1 ] )
                %lightangle(-45,30)
                hold on
                plot3(storedpath(:,1),storedpath(:,2),storedpath(:,3))
                
                plot3(storedpath(end-1:end,1),storedpath(end-1:end,2),storedpath(end-1:end,3),'r')
                xlim([-cellradius-0.1 ,cellradius+0.1 ]); ylim([-cellradius-0.1 ,cellradius + 0.1 ]); zlim([-cellradius-0.1 ,cellradius + 0.1 ]);
                xlabel('x');ylabel('y');zlabel('z');
                axis equal
                drawnow
                set(gcf,'color','w');
                frame=getframe(gcf);
                writeVideo(videowrite,frame);
                if numbsteps == 1000, close(videowrite); end
                drawnow
                %pause(0.5)
            end
        end
        
        allangles = []; %store all angles here
        diffvects  = storedpath(2:end,:) - storedpath(1:end-1,:);
        
        for j = 2:size(diffvects,1)
            
            diffang = atan2(norm(cross(diffvects(j,:),diffvects(j-1,:))), dot(diffvects(j,:),diffvects(j-1,:))); %find difference in angle
            allangles = [allangles, diffang];
            
        end
        
        %histogram(allangles)
        
        if strcmp(filename,'SampleFromUniformAngle')
            save(['insilico_' , filename ],'allangles','diffvects','storedpath','vectr');
        else
            save(['insilico_' , filename ],'allangles','diffvects','vectr','storedpath','allx','ally','diffangs');
        end
    end
end

%% Generate graphs of persistence

if plot_persistence == 1
    
    %% Compare persistence between in vivo data and in silico stuff
    
    %generate persistence graphs from in vivo data
    
    allangles = []; allrs = []; %store all persistence angles and step sizes here
    
    for cellfilenames = {'AllTracks_Data from 20160206_PDGGTexp3_18172_18274_Spa2mCh0nM_dly18172.mat' , 'AllTracks_Data from 20160511_PDGGTexp11_18172_211090nM_dly18172.mat' }
        %for cellfilenames = {'AllTracks_Data from Allie 140330_dose_msds_1817218172_0nM.mat'}
        %for cellfilenames =  {'AllTracks_Data from 20160206_PDGGTexp3_18172_18274_Spa2mCh0nM_dly18172.mat' , 'AllTracks_Data from 20160511_PDGGTexp11_18172_211090nM_dly18172.mat','AllTracks_Data from Allie 140330_dose_msds_1817218172_0nM.mat' }
        %for cellfilenames =  {'AllTracks_Data from Allie 150921_latAlata.mat'}
        
        filename = char(cellfilenames);
        load(filename);
        
        for iindex = 1:size(alltracks,2) %alltracks is a cell matrix with all the tracks
            
            diffvects  = alltracks{iindex}(2:end,2:4) - alltracks{iindex}(1:end-1,2:4); % movement vectors of all the patches
            
            for j = 2:size(diffvects,1)
                
                diffang = atan2(norm(cross(diffvects(j,:),diffvects(j-1,:))), dot(diffvects(j,:),diffvects(j-1,:))); %find difference in angle
                allangles = [allangles, diffang];
                
            end
            
            allrs = [allrs ; sqrt(diffvects(:,1).^2 + diffvects(:,2).^2 + diffvects(:,3).^2) ];
            
        end
        
    end
    %save('lata.mat','allrs','allangles')
    
    % plot persistence values and step lengths from experiments
    
    figure
    
    subplot(1,2,1)
    histogram(allrs,'BinWidth',0.1)
    axis square
    xlim([0 4])
    
    subplot(1,2,2)
    histogram(allrs)
    axis square
    xlim([0 4])
    set(gca,'YScale','log') %plot on log log scale
    set(gca,'XScale','log')
    
    figure %plot data generated from in vivo metrics
    subplot(1,2,1)
    hold on
    cdfplot(allangles); xlim([0,pi]);
    
    subplot(1,2,2)
    hold on
    histogram(allrs,'BinWidth',0.1,'Normalization','probability','DisplayStyle','stairs')
    axis square
    xlim([0 4])
    
    
    % plot persistence values and step length distribution from in silico
    % model
    
    for filename = {'insilico_20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024.mat'}
        load(char(filename));
        
        subplot(1,2,1) % plot persistence
        hold on
        cdfplot(allangles);
        
        axis square
        xlabel('Theta (rad)')
        
        subplot(1,2,2) % plot distance traveled per unit length
        hold on
        cdfplot(vectr(:))
        %histogram(vectr(:),'BinWidth',0.1,'Normalization','probability','DisplayStyle','stairs')
        axis square
        xlim([0 4])
        %legend('exp','0','1','2','3','5')
        
        %set(gca,'YScale','log') %plot on log log scale
        %set(gca,'XScale','log')
        
        set(gcf ,'renderer' ,'painters') %do this or matlab will fucking fail to generate a vector graphic
        
    end
    
    
    %% plot changes in persistence across different uniform pheromone concentrations
    figure
    
    allnames = {{'AllTracks_Data from 20160206_PDGGTexp3_18172_18274_Spa2mCh0nM_dly18172.mat' , 'AllTracks_Data from 20160511_PDGGTexp11_18172_211090nM_dly18172.mat'},...
        {'AllTracks_Data from 20160206_PDGGTexp3_18172_18274_Spa2mCh20nM_dly18172.mat'},...
        {'AllTracks_Data from 20160206_PDGGTexp3_18172_18274_Spa2mCh40nM_dly18172.mat'},...
        {'AllTracks_Data from 20160206_PDGGTexp3_18172_18274_Spa2mCh60nM_dly18172.mat'},...
        {'AllTracks_Data from 20160206_PDGGTexp3_18172_18274_Spa2mCh80nM_dly18172.mat'}};
    
    for pherconc = 1:5
        
        allangles = []; allrs = []; %store all persistence angles and step sizes here
        
        for withinphercon = 1:size(allnames{pherconc},2)
            
            filename = char(allnames{pherconc}{withinphercon});
            load(filename);
            
            for iindex = 1:size(alltracks,2)
                
                diffvects  = alltracks{iindex}(2:end,2:4) - alltracks{iindex}(1:end-1,2:4);
                
                for j = 2:size(diffvects,1)
                    
                    diffang = atan2(norm(cross(diffvects(j,:),diffvects(j-1,:))), dot(diffvects(j,:),diffvects(j-1,:))); %find difference in angle
                    allangles = [allangles, diffang];
                    
                end
                
                allrs = [allrs ; sqrt(diffvects(:,1).^2 + diffvects(:,2).^2 + diffvects(:,3).^2) ];
                
            end
            
        end
        
        subplot(1,2,1)
        hold on
        cdfplot(allangles);
        xlabel('Theta (rad)')
        axis square
        
        subplot(1,2,2)
        hold on
        cdfplot(allrs);
        xlabel('Distance (micron)')
        axis square
        
    end
    
    legend('0','20','40','60','80');
    set(gcf ,'renderer' ,'painters') %do this or matlab will fucking fail to generate a vector graphic
    
end

%% find autocorrelation for persistent behavior

if persistenceautocorrelation == 1 %see what autocoreelation graph looks like for persistent movement
    
    amp = 1;
    sigmacorr = 0.9;
    uni = 0.2;
    
    
    simnum = 1000000;
    stepnum = 10;
    
    
    allpersdir = zeros(simnum,stepnum);
    
    for ii = 1:simnum
        
        for jj = 2:stepnum
            
            newang = allpersdir(ii,jj-1) + Function_sample_normal_uniform( amp , sigmacorr , uni );
            
            allpersdir(ii,jj) = mod(newang - pi , 2*pi) - pi; %constrain within +/-pi
            
        end
        
    end
    
    allrandangles = rand(simnum,stepnum)*2*pi - pi;
    allrandangles(:,1) = 0;
    
    %see how signal decay looks
    
    allperssig = (pi/2 - mean(abs(allpersdir),1))/(pi/2);
    allrandsig = (pi/2 - mean(abs(allrandangles),1))/(pi/2);
    
    hold on
    timeaxis = [0:stepnum-1]*1.5;
    plot(timeaxis,allperssig) %half life tau = 1.0383
    plot(timeaxis,allrandsig)
    plot(timeaxis,allperssig,'b.')
    plot(timeaxis,allrandsig,'r.')
    
    
    
    axis square
    
    %{
    corrdir = xcorr(allpersdir,allpersdir);
    plot(corrdir)
    hold on
    corrrand = xcorr(allrandangles,allrandangles);
    plot(corrrand)
    %autocorr(allpersdir,'NumLags',100)
    %}
    
end

if persistenceaudocorrelation_experimental_data == 1
    
    load('Invivo.mat');
    
end

%% plot step size and persistence for different situations, like planar vs spherical, different cell radii
%
%      Persistence Angle CDF           Step length histogram
%
%                                      x
%      x                xxx            x      xx
%      x           xxxxxx              x     xxx
%      x        xxx                    x    xx xx
%      x      xxx                      x    x   x
%      x     xx                        x   xx   xx
%      x    xx                         x   x     xx
%      x   xx                          x   x      xxxx
%      x  xx                           x xxx         xxx
%      x  x                            xxx             xxxx
%      x x                             xx                  xx
%      xxxxxxxxxxxxxxxxxxxxx           xxxxxxxxxxxxxxxxxxxxxxxxxx
%              angles                          distances

if plot_persistence_2 == 1
    
    % plot persistence for experimental stuff
    
    threshholdstep = 1.5; %step lengths below which you want to compare with the model
    cdfunifarea = 1.5493; %area of CDF of angles for uniform sampling of angles -- thus, with no persistence
    showtext = 1; %see if you want to show text labels in figure
    plotcount = 1; %keep track of how many points you plotted
    
    load('Invivo.mat');
    allrs = allrs(allrs<threshholdstep);
    figure
    subplot(2,3,1) %plot angle cdfs
    
    hold on
    
    angprop = cdfplot(allangles);
    cdfangarea = trapz(angprop.XData(2:(end-1)),angprop.YData(2:(end-1))); % Calculate total trapezoidal area
    cdfangdiff = cdfangarea - cdfunifarea; %caclulate deviation from uniform sampling of angles while moving on a sphere -- that is, difference in areas of CDFs
    
    allanglesexp = allangles; %this is so you can keep the experimental angles for kstest later
    
    xlim([0,pi]); %plot CDF of persistence angles
    
    subplot(2,3,2) %plot step length distributions
    
    hold on
    [histvalues,histedges] = histcounts(allrs,'BinWidth',0.1,'Normalization','probability');
    
    area(histedges(1:(end-1)),histvalues(1:end))
    axis square
    
    xlabel('Step length (\mum)')
    
    disp( ['Mean exp:' ,    num2str( mean(allrs(allrs<threshholdstep)) ) ]) %find mean of just rs that are less than 2um
    disp([ 'Std exp: '  , num2str( std(allrs(allrs<threshholdstep)) ) ])
    
    subplot(2,3,3) %plot mean vs std of step lengths
    
    hold on
    plot( std(allrs(allrs<threshholdstep)) , mean(allrs(allrs<threshholdstep))  ,'.' ,'MarkerSize', 20) ;
    
    subplot(2,3,4) %d-stats of step-length vs turning angles
    
    hold on
    plot( 0  , 0 , '.' , 'MarkerSize', 20);
    
    
    subplot(2,3,5) % 3d plot of mean, std, and persistence
    
    hold on
    plot( cdfangdiff , mean(allrs(allrs<threshholdstep)) , '.' ,'MarkerSize' , 20  )
    
    xlabel('Pers'); ylabel('Mean');
    
    subplot(2,3,6) % std vs pers
    
    hold on
    plot( std(allrs(allrs<threshholdstep)) ,   cdfangdiff , '.' ,'MarkerSize' , 20  )
    
    xlabel('Std'); ylabel('Pers');
    
    % straight up comparison between base model and experiments
    
    %{
        
         for filename = {'SampleFromUniformAngle','20190113_SETMAWP_Treadmill_0_0'
            }
        
             legendstuff = {'Exp','Unif','Base'};
             
             meanlow = 0.2; meanhigh = 0.5;
             perslow = 0.15; pershigh = 0.5;
             stdlow = 0.1; stdhigh = 0.3;
        
    %}
    
    % All three -- focused cables, temporal variation, and local GAP
    
    %{%
    
    for filename = {'SampleFromUniformAngle','20190113_SETMAWP_Treadmill_0_0',...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_5_actinvar_024_LocalGAP_thalf_10_GAP_0_3_clusteredvesicle_LtH_100_HtL_20_LR_50'
            }
        
        legendstuff = {'Exp','Unif','Base','Full Model'};
    
    
    %}
    
    
     % Exocytosis and endocytosis made uniform
    
     %{
     
     for filename = {'SampleFromUniformAngle','20190113_SETMAWP_Treadmill_0_0',...
             '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_0_actinvar_024_unif_endocytosis',...
             '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_0_actinvar_024_unif_exocytosis'}
         
         legendstuff = {'Exp','Unif','Base','Endo','Exo'};
        
        
    
    % Bem-Snc simulation
    
     %}
     
     %{
     
     for filename = {'SampleFromUniformAngle','20190113_SETMAWP_Treadmill_0_0',...
             '20190113_SETMAWP_Treadmill_0_0_BemSnc1_k3v_0_30',...
             '20190113_G_SETMAWP_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_clusteredvesicle_LtH_100_HtL_20_LR_100_k3_06645'}
         
         legendstuff = {'Exp','Unif','Base','Bem-Snc','High GEF'};
         
         
         %}
         
         %increasing vesicular GAP activity
         
         %{
        
        for filename = {'SampleFromUniformAngle',...
                '20190113_SETMAWP_Treadmill_0_0',...
                '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_0_actinvar_024_LocalGAP_thalf_5_GAP_0_5_clusteredvesicle_LtH_60_HtL_20_LR_100',...
                '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_5_clusteredvesicle_LtH_60_HtL_20_LR_100',...
                '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_0_actinvar_024_LocalGAP_thalf_5_GAP_1_clusteredvesicle_LtH_60_HtL_20_LR_100'
                }
            
            legendstuff = {'Exp','Unif','Base','5s 0.5','10s 0.5','5s 1'};
            
            %}
        
        %changing focusedness of actin cables
        
        %{ 
        
        for filename = {'SampleFromUniformAngle',...
                '20190113_SETMAWP_Treadmill_0_0',...
                '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024',...
                '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_3_actinvar_024',...
                '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_5_actinvar_024'...
                }
            legendstuff = {'Exp','Unif','kG=0','kG=1','kG=3','kG=5'};
            
            %}
        
        
        %changing diffusion constant
        
        %{
    
    
    for filename = {'SampleFromUniformAngle','20190113_SETMAWP_Treadmill_0_0','20190113_SETMAWP_Treadmill_0_0_Diff_0005', '20190113_SETMAWP_Treadmill_0_0_Diff_01'...
            }
        legendstuff = {'Exp','Unif','0.0045','0.0005','0.01'};
        %}
        
        %{
    
    %changing length scale of clustering
    
    for filename = {'SampleFromUniformAngle','20190113_SETMAWP_Treadmill_0_0',...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_clusteredvesicle_LtH_25_HtL_5_LR_60',...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_clusteredvesicle_LtH_50_HtL_10_LR_60',...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_clusteredvesicle_LtH_100_HtL_20_LR_60',...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_clusteredvesicle_LtH_200_HtL_40_LR_60'
            }
        legendstuff = {'Exp','Unif','Base','5:25','10:50','20:100','40:200'};
        %}
        

        %testing MWE vs regular base model
        
        %{
        
        for filename = {'SampleFromUniformAngle','20190113_SETMAWP_Treadmill_0_0',...
        '20190113_G_SETMAWP_k_gauss_1_5_actinvar_024_LocalGAP_thalf_10_GAP_0_clusteredvesicle_LtH_100_HtL_20_LR_100',...
            '20190113_G_SETMAWP_k_gauss_1_5_actinvar_024_LocalGAP_thalf_10_GAP_0_3_clusteredvesicle_LtH_100_HtL_20_LR_50',...
            '20190113_G_SETMAWP_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_clusteredvesicle_LtH_100_HtL_20_LR_50',...
            '20190113_G_SETMAWP_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_3_clusteredvesicle_LtH_100_HtL_20_LR_100',...
            '20190113_G_SETMAWP_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_clusteredvesicle_LtH_100_HtL_20_LR_100'}
        
        legendstuff = {'Exp','Unif','Base','kG=1.5','Full','Burst=50%','GAP=0.3','MWE'};
        
        
        %}
        
        %changing basal vesicle delivery rate
        
        %{
    
    for filename = {'SampleFromUniformAngle','20190113_SETMAWP_Treadmill_0_0',...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_clusteredvesicle_LtH_100_HtL_20_LR_80',...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_clusteredvesicle_LtH_100_HtL_20_LR_60',...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_clusteredvesicle_LtH_100_HtL_20_LR_40'
            }
        legendstuff = {'Exp','Unif','Base','80','60','40'};
    
        %}
        
        
        %changing vesicular GAP effect
        
        %{
        
         for filename = {'SampleFromUniformAngle','20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_clusteredvesicle_LtH_60_HtL_20_LR_60',...
                 '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_LocalGAP_thalf_10_GAP_0_5_clusteredvesicle_LtH_60_HtL_20_LR_60',...
                 '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_LocalGAP_thalf_10_GAP_1_clusteredvesicle_LtH_60_HtL_20_LR_60'}
        
             legendstuff = {'Exp','Unif','0','0.5','1'};
    
    meanlow = 0.2; meanhigh = 0.4;
    perslow = 0; pershigh = 1;
    stdlow = 0.1; stdhigh = 0.4;
    
        %}
        
        
        
        % Changing overall GAP acitvity and diffusion constant
        
        %{
    for filename = {'SampleFromUniformAngle','20190113_SETMAWP_Treadmill_0_0',...
            '20190113_SETMAWP_Treadmill_0_0_k2b_0_7',...
            '20190113_SETMAWP_Treadmill_0_0_k2b_0_9',...
            '20190113_SETMAWP_Treadmill_0_0_Diff_0005',...
            '20190113_SETMAWP_Treadmill_0_0_Diff_0015',...
            '20190113_SETMAWP_Treadmill_0_0_Diff_01'
            }
        
             legendstuff = {'Exp','Unif','Base','0.7','0.9','0.0005','0.0015','0.01'};
             
             
        %}
        
        %Changing vesicle delivery rate
        
        %{
        
        for filename = {'SampleFromUniformAngle','20190113_SETMAWP_Treadmill_0_0',...
                '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_periodOut_03',...
                '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_periodOut_01'}
            
               legendstuff = {'Exp','Unif','Base','2X','6X'};   
        
        %}
        
        meanlow = 0.2; meanhigh = 0.45;
        perslow = 0; pershigh = 1;
        stdlow = 0.1; stdhigh = 0.3;
        
        plotcount = plotcount + 1;
        
        load(['insilico_',char(filename),'.mat']);
        
        cellradius = 2.5;
        
        subplot(2,3,1) % plot persistence
        
        hold on
        angprop = cdfplot(allangles);
        cdfangarea = trapz(angprop.XData(2:(end-1)),angprop.YData(2:(end-1))); % Calculate total trapezoidal area
        cdfangdiff = cdfangarea - cdfunifarea; %caclulate deviation from uniform sampling of angles while moving on a sphere -- that is, difference in areas of CDFs
        
        if strcmp('SampleFromUniformAngle',char(filename)) %just to find out what the area is for the uniform case
            cdfunifarea = cdfangarea; % found that cdfunifarea is  1.5493
        end
        
        axis square
        grid off
        xlabel('Theta (rad)')
        
        subplot(2,3,2) % plot distance traveled per unit length
        
        hold on
        sphd = cellradius*sqrt(2*(1-cos(vectr(:)/cellradius)));
        
        [histvalues,histedges] = histcounts(sphd,'BinWidth',0.1,'Normalization','probability');
        
        plot(histedges(1:(end-1)),histvalues(1:end))
        axis square
        xlim([0 threshholdstep])
        xlabel('Step length (\mum)')
        
        %Angle statistics
        disp('----------------')
        disp([ filename])
        
        [ hypoth , pvalu , kstatsang  ] = kstest2(allanglesexp,allangles);
        disp( ['Angp = ' , num2str(pvalu) , ' AngleDStat = ' , num2str(kstatsang) , ' Hyp = ' , num2str(hypoth)  ] );
        
        % Step length statistics
        
        [ hypoth , pvalu , kstatslen  ] = kstest2(allrs(allrs<threshholdstep),sphd);
        disp( ['Steppvalu = ' , num2str(pvalu) , ' StepDStat = ' , num2str(kstatslen) , ' Hyp = ' , num2str(hypoth)  ] );
        
        disp( ['Meanstep:' ,    num2str( mean(sphd(:)) ) ])
        disp([ 'Stdstep: '  , num2str( std(sphd(:)) ) ])
        
        subplot(2,3,3) %plot mean vs std of step lengths
        
        hold on
        plot( std(sphd(:)) , mean(sphd(:)) , '.' ,'MarkerSize', 20);
        xlabel('Std'); ylabel('Mean')
        axis square
        
        if showtext == 1 %label point
            text(std(sphd(:)) , mean(sphd(:)) , char(legendstuff(plotcount)) )
            
        end
        
        
        xlim([  stdlow stdhigh ]) %std limits
        ylim([meanlow meanhigh]) %Mean limits
        
        
        subplot(2,3,4) %d-stats of step-length vs turning angles
        
        hold on
        plot( kstatslen  , kstatsang , '.' , 'MarkerSize', 20);
        
        xlabel('Dstatlen'); ylabel('Dstatang')
        axis square
        
        
        subplot(2,3,5) % mean vs pers
        
        hold on
        plot(  cdfangdiff , mean(sphd(:)) ,  '.' ,'MarkerSize' , 20  )
        
        xlabel('Pers'); ylabel('Mean');
       
        xlim([ perslow pershigh ]) %persistence limits
        ylim([meanlow meanhigh]) %Mean limits
        
        if showtext == 1 %label point
            text( cdfangdiff , mean(sphd(:))  , char(legendstuff(plotcount)) )
            
        end
        
        axis square
        
        subplot(2,3,6) % std vs pers
        
        hold on
        plot( std(sphd(:)) ,   cdfangdiff , '.' ,'MarkerSize' , 20  )
        
        xlabel('Std'); ylabel('Pers');
        
        xlim([ stdlow stdhigh ]) %std limits
        ylim([ perslow pershigh ]) %persistence limits
        
        if showtext == 1 %label point
            text(  std(sphd(:)) ,   cdfangdiff , char(legendstuff(plotcount)) )
            
        end
        
        axis square
        
    end
    
    legend(legendstuff)
    
    set(gcf, 'Position', [0 0 740 430])
    
    
end

%% plot step size and persistence for different situations, like planar vs spherical, different cell radii
%
%      Persistence Angle CDF           Step length histogram
%
%                                      x
%      x                xxx            x      xx
%      x           xxxxxx              x     xxx
%      x        xxx                    x    xx xx
%      x      xxx                      x    x   x
%      x     xx                        x   xx   xx
%      x    xx                         x   x     xx
%      x   xx                          x   x      xxxx
%      x  xx                           x xxx         xxx
%      x  x                            xxx             xxxx
%      x x                             xx                  xx
%      xxxxxxxxxxxxxxxxxxxxx           xxxxxxxxxxxxxxxxxxxxxxxxxx
%              angles                          distances

if plot_persistence_3 == 1
    
    load('AllModelsNames.mat') %load names and plotting color for each model
    
    % plot persistence for experimental stuff
    
    alldistmeanstd = [];
    alldistangk = [];
    
    threshholdstep = 1.5; %step lengths below which you want to compare with the model
    
    load('Invivo.mat');
    allrs = allrs(allrs<threshholdstep);
    figure
    subplot(2,2,1) %plot angle cdfs
    
    hold on
    allanglesexp = allangles; %this is so you can keep the experimental angles for kstest later
    cdfplot(allanglesexp); xlim([0,pi]); %plot CDF of persistence angles
    
    subplot(2,2,2) %plot step length distributions
    
    hold on
    [histvalues,histedges] = histcounts(allrs,'BinWidth',0.1,'Normalization','probability');
    
    plot(histedges(2:end),histvalues(1:end))
    axis square
    xlim([0 2])
    xlabel('Step length (\mum)')
    
    disp( ['Mean exp:' ,    num2str( mean(allrs(allrs<threshholdstep)) ) ]) %find mean of just rs that are less than 2um
    disp([ 'Std exp: '  , num2str( std(allrs(allrs<threshholdstep)) ) ])
    
    subplot(2,2,3) %plot mean vs std of step lengths
    
    hold on
    plot( mean(allrs(allrs<threshholdstep)) ,  std(allrs(allrs<threshholdstep)) ,'.' ,'MarkerSize', 20 , 'MarkerEdgeColor', AllModels{1,1}) ;
    
    alldistmeanstd = [ alldistmeanstd ; mean(allrs(allrs<threshholdstep)) ,  std(allrs(allrs<threshholdstep))  ];
    
    
    subplot(2,2,4) %d-stats of step-length vs turning angles
    
    hold on
    plot( 0  , 0 , '.' , 'MarkerSize', 20 , 'MarkerEdgeColor', AllModels{1,1});
    
    alldistangk = [alldistangk ; 0 , 0];
    
    % plot persistence for in silico stuff
    
    dotcolors = [ 0 0 0;... %exp
        1 1 1;...  unif
        1 0 0;...  base
        1 0 0;...  diff 1
        1 0 0;...  diff 2
        0 0 1;...  kG = 1
        0 0 1;...  kG = 2
        0 0 1;...  kG = 3
        0 0 1;...  kG = 5
        0 1 1;...  kG = 1; GAP = 090
        1 0 1;...  kg = 1; LtH_100_HtL_20
        1 0 1;
        1 0.6 0;...6X vesicles delivery
        1 0.6 0;...
        1 0.6 0;...
        1 0.6 1;...GAP + LtH_60_HtL_20 + kG = 1
        ];
    
    
    for ii = 2:size(AllModels,1)
        
        filename = AllModels{ii,2};
        
        %{
        filename = {'SampleFromUniformAngle' ...  unif
            '20190113_SETMAWP_Treadmill_0_0' ...  base
            '20190113_SETMAWP_Treadmill_0_0_Diff_0025' ... diff 1
            '20190113_SETMAWP_Treadmill_0_0_Diff_0060' ... diff 2
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024'... kg = 1
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_2_actinvar_024'...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_3_actinvar_024'...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_5_actinvar_024'... kg = 5
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_GAP_090'... kg = 1; GAP = 090
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_clusteredvesicle_LtH_100_HtL_20_LR_60'... kg = 1; LtH_100_HtL_20
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_5_actinvar_024_clusteredvesicle_LtH_60_HtL_20_LR_60'...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_periodOut_01'... 6X vesicles delivery
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_periodOut_03'...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_periodOut_08'...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_LocalGAP_thalf_10_GAP_1_clusteredvesicle_LtH_60_HtL_20_LR_60'...GAP + LtH_60_HtL_20 + kG = 1
            }
        %}
        
        
        % for filename = {'SampleFromUniformAngle','20190113_SETMAWP_Treadmill_0_0','20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024',...
        %         '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_clusteredvesicle_LtH_100_HtL_20_LR_60'}
        
        load(['insilico_',char(filename),'.mat']);
        
        cellradius = 2.5;
        
        subplot(2,2,1) % plot persistence
        
        hold on
        cdfplot(allangles);
        
        axis square
        xlabel('Theta (rad)')
        
        subplot(2,2,2) % plot distance traveled per unit length
        
        hold on
        sphd = cellradius*sqrt(2*(1-cos(vectr(:)/cellradius)));
        
        [histvalues,histedges] = histcounts(sphd,'BinWidth',0.1,'Normalization','probability');
        
        plot(histedges(2:end),histvalues(1:end),'Color',AllModels{ii,1})
        
        axis square
        xlim([0 2])
        xlabel('Step length (\mum)')
        
        
        %Angle statistics
        disp('----------------')
        disp([ filename])
        
        [ hypoth , pvalu , kstatsang  ] = kstest2(allanglesexp,allangles);
        disp( ['Angp = ' , num2str(pvalu) , ' AngleDStat = ' , num2str(kstatsang) , ' Hyp = ' , num2str(hypoth)  ] );
        
        % Step length statistics
        
        [ hypoth , pvalu , kstatslen  ] = kstest2(allrs(allrs<threshholdstep),sphd);
        disp( ['Steppvalu = ' , num2str(pvalu) , ' StepDStat = ' , num2str(kstatslen) , ' Hyp = ' , num2str(hypoth)  ] );
        
        disp( ['Meanstep:' ,    num2str( mean(sphd(:)) ) ])
        disp([ 'Stdstep: '  , num2str( std(sphd(:)) ) ])
        
        
        subplot(2,2,3) %plot mean vs std of step lengths
        
        hold on
        plot( mean(sphd(:)) ,  std(sphd(:)) , '.' ,'MarkerSize', 20 , 'MarkerEdgeColor', AllModels{ii,1});
        xlabel('Mean'); ylabel('Std')
        axis square
        
        alldistmeanstd = [ alldistmeanstd ;  mean(sphd(:)) ,  std(sphd(:)) ];
        
        
        subplot(2,2,4) %d-stats of step-length vs turning angles
        
        hold on
        plot( kstatslen  , kstatsang , '.' , 'MarkerSize', 20, 'MarkerEdgeColor', AllModels{ii,1});
        
        xlabel('Dstatlen'); ylabel('Dstatang')
        axis square
        
        
        alldistangk = [ alldistangk ; kstatslen  , kstatsang  ];
        
    end
    
    
    %{
    figure
    
   subplot(1,2,1)
   scatter( alldistmeanstd(:,1) , alldistmeanstd(:,2), 20,dotcolors, 'filled')
   xlabel('Step length mean')
   ylabel('Step length std')
   axis square
   
   subplot(1,2,2)
   scatter( alldistangk(:,1) , alldistangk(:,2), 20,dotcolors, 'filled')
   axis square
   xlabel('kstat Step'); ylabel('kstat Angle')
    %}
end

%% generate graphs of msd or diffusion constants for the simulations and experiments

if plot_diffusion_constants == 1
    
    
    %% first do this for in vivo experiments
    
    timeaxis = [1:20]*1.5; %time in minutes
    sumsds = zeros(1,20); %takes sum of all square displacements
    sdno = zeros(1,20); %how many square displacements are input into sumsds
    
    for cellfilenames = {'AllTracks_Data from 20160206_PDGGTexp3_18172_18274_Spa2mCh0nM_dly18172.mat' , 'AllTracks_Data from 20160511_PDGGTexp11_18172_211090nM_dly18172.mat' }
        %for cellfilenames = {'AllTracks_Data from Allie 140330_dose_msds_1817218172_0nM.mat'}
        %for cellfilenames =  {'AllTracks_Data from 20160206_PDGGTexp3_18172_18274_Spa2mCh0nM_dly18172.mat' , 'AllTracks_Data from 20160511_PDGGTexp11_18172_211090nM_dly18172.mat','AllTracks_Data from Allie 140330_dose_msds_1817218172_0nM.mat' }
        
        filename = char(cellfilenames);
        load(filename)
        
        for trackno = 1:size(alltracks,2) %pick track number
            
            for windowsize = 1:(size(alltracks{trackno},1)-1) % window size to calculate msd
                
                for startstep = 1:(size(alltracks{trackno},1) - windowsize)
                    
                    init = alltracks{trackno}(startstep,2:4); fina = alltracks{trackno}(startstep+windowsize,2:4); %initial and final positions of the tracks
                    sumsds(windowsize) = sumsds(windowsize) + norm(fina - init)^2; %find square displacement
                    sdno(windowsize) = sdno(windowsize)+1;
                    
                end
            end
        end
    end
    
    msds = sumsds./sdno;
    hold on
    plot([0,timeaxis],[0,msds],'b');
    
    %% do the same for simulation of movement of in silico path on sphere
    
    msds = []; %where msds are stored
    
    for cellfilenames = {'insilico_gauss_1_actinvar_024.mat' }
        filename = char(cellfilenames);
        load(filename)
        
        for windowsize = 1:20 % window size to calculate msd
            
            sds = storedpath(1:windowsize:end,:);
            sds = sum((sds(2:end,:) - sds(1:(end-1),:)).^2,2); %square distance traveled
            msds(windowsize) = mean(sds(:));
            
        end
        
        hold on
        plot([0,timeaxis],[0,msds]);
        
    end
    
    %% do the same with patch diffusing along flat surface
    
    set(gca,'ColorOrderIndex',1) %reset color order
    
    msds = []; %where msds are stored
    
    for cellfilenames = {'insilico_gauss_1_actinvar_024.mat' }
        filename = char(cellfilenames);
        load(filename)
        
        for windowsize = 1:20 % window size to calculate msd
            
            sdsx = allx(:,1:windowsize:end); sdsy = ally(:,1:windowsize:end);
            sds = (sdsx(:,2:end) - sdsx(:,1:(end-1))).^2 + (sdsy(:,2:end) - sdsy(:,1:(end-1))).^2 ;
            msds(windowsize) = mean(sds(:));
            
        end
        
        hold on
        plot([0,timeaxis],[0,msds],'--');
        
    end
    
    
    %{
    %% include controls for debugging purposes
    
    storedpath = zeros(10000,2);
    prevang = 2*pi*rand();
    for i = 2:size(storedpath,1)
        
        %storedpath(i,:) = storedpath(i-1,:) + (rand(1,2)-0.5); %diffusive motion
        
        
    end
    
    for windowsize = 1:20 % window size to calculate msd
        
        sds = storedpath(1:windowsize:end,:);
        sds = sum((sds(2:end,:) - sds(1:(end-1),:)).^2,2); %square distance traveled
        msds(windowsize) = mean(sds(:));
        
    end
    
    hold on
    plot([0,timeaxis],[0,msds],'b');
    %}
    
    %% plot properties
    
    legend('In vivo','0','1','5','Flat0','Flat1','Flat5')
    xlim([0,10.5]);
    xlabel('Minutes')
    ylabel('Msd')
    
    axis square
    
end


%% Simulate patch movment on spheres that are near each other
%
%         Adjacent wandering cells                                    Time
%                                                                              XXX
%                                                                             XXX            +
%                                                              +           XXX               |
%      XXXXXXXXXXXXX           XXXXXXXXX                       |          XX                 |        X
%    XXX           XX       XXX        XXX                     |         X                   |
%  X X              XX     XX            XXX                   |        X                    |
% XX                 XX    X               X                   |      XX        XXX          |
% X             XXXX   X  XX               XXX                 |     XX      XXXX            |
% X             XXXXX  X  X      XX  XXXX    X                 |    XX     XXX               |                 X
% X   XXX      XXXXX  X   X    XXXXX X   X   XX                |   XX   XXX                  |   X
% XXXX   XXXXXX       X   XX    XXX XX   X    X                |  XX  XX         XXXXX       |
%  XX                 X    X           XXX   XX                | X  XX      XXXXX            |
%  XX               XXX    XX        XX     XX                 |X X X   XXXX                 |
%    XX            XX       XXX      X    XX                   |XXXXX                        |
%     XXXXXXXXXXXXX           XXXXXXXXXXXX                     +-----------------------+     +----------------------+
%                                                                         time                        angles

if simulate_search_efficiency == 1
    
    simulationtype = 3;
    % 1. simulate patch movement by sampling/modeling r and theta from offset Gaussian
    % 2. simulate patch movement by modeling theta from experiment but r from a gaussian that was generated by fitting experimental msd
    % 3. simulate patch movement by sampling r and theta from the large PDE model
    
    if simulationtype == 1
        
        myparameters = [0.1 0.3 0.5 0.9 1.5 2 3];        
        
    elseif simulationtype == 2
        
        myparameters = [0];
        
    elseif simulationtype == 3
        
         myparameters = {'insilico_20190113_G_SETMAWP_k_gauss_1_5_actinvar_024_LocalGAP_thalf_10_GAP_0_clusteredvesicle_LtH_100_HtL_20_LR_100',...
            'insilico_20190113_G_SETMAWP_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_clusteredvesicle_LtH_100_HtL_20_LR_50',...
            'insilico_20190113_G_SETMAWP_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_3_clusteredvesicle_LtH_100_HtL_20_LR_100',...
            'insilico_20190113_G_SETMAWP_k_gauss_0_actinvar_024_LocalGAP_thalf_10_GAP_0_clusteredvesicle_LtH_100_HtL_20_LR_100',...
            'insilico_20190113_G_SETMAWP_k_gauss_1_5_actinvar_024_LocalGAP_thalf_10_GAP_0_3_clusteredvesicle_LtH_100_HtL_20_LR_50'};
    end
    
    
    for parameter =  myparameters %sweep through some parameter -- you can change what the parameter is
        
       
        cellradius = 2.5; %microns
        startpositions = [ 1 , -1 ]; %parameters that determine two positions where you want the patches to start
        numbsteps = 100; %number of simulations steps you want this to take
        numbsims = 5000; %number of cells you want to simulate
        timeinterval = 1.5; %minutes
        
        sphere_paths = cell(numbsims,numel(startpositions));
        
        if simulationtype == 1 || simulationtype == 2
            
            load('InVivo.mat');
            convertedrs = real(cellradius*acos(1 - allrs.^2/(2*cellradius^2))); % since these are all straight lines, convert to how far a point would have traveled on a sphere of some radius
            
            
            if simulationtype == 1 %establish values for amp, sigma, and uni
                %amp = 1;
                %sigma = 0.9;
                uni = 0.2;
                
                normfunctionexp = @(x) 1*exp( -(x).^2/(2*0.9^2) ); %use this to find area of gaussian part of 'offset gaussian' that fits experimental data.
                %Recall that  amp = 1, sigma = 0.9, and uni = 0.2 seems to match turning angle distribution seen in experiments.
                angareaexp = integral(normfunctionexp,-pi,pi);
                
                normfunction = @(x) 1*exp( -(x).^2/(2*parameter^2) ); %use this to calibrate 'amp' so that area is the same
                angarea = integral(normfunction,-pi,pi);
                
                amp = angareaexp/angarea; %new amp value
                
            end
            
        elseif simulationtype == 3
            
            load(char(parameter));
            
        end
        
        %convertedrs = allrs;
        if 1 %open flag if you want to simulate stuff
            
            for sim = 1:numbsims
                
                for startat = 1:numel(startpositions) %start with patch on back half of each spherical cell
                    
                    storedpath = zeros(numbsteps,3); % initialize martix where you will store patch path
                    
                    %generate a point on the back half of spherical cell
                    
                    ele = pi*rand() + startpositions(startat)*pi/2;
                    az = asin(2*rand()-1);
                    [Xc,Yc,Zc] = sph2cart(ele,az,cellradius);
                    
                    storedpath(1,:) = [ Xc,Yc,Zc ];
                    
                    
                    for iindex = 2:numbsteps
                       
                        if simulationtype == 1
                            
                            %simulate experimental cell's movement with both r and theta
                            ang = Function_sample_normal_uniform( amp , parameter , uni ); % amp = 1, sigma = 0.9, and uni = 0.2 seems to match turning angle distribution seen in experiments
                            mydist = randsample(convertedrs,1);
                            
                        elseif simulationtype == 2
                            
                            %simulate same thing but using theta from experiments and r from a gaussian that was chosen by fitting msd to experimental data
                            ang = rand()*2*pi-pi; %sample with no persistence
                            mydist = abs(normrnd(0,sqrt(4*0.0971*timeinterval))); %0.0971 um^2/min is the diffusion constant extracted from simulated persistent movement from in vivo data
                            
                        elseif simulationtype == 3
                            
                            %simulate by sampling r and theta from the PDE
                            ang = randsample(diffangs(:),1);
                            mydist = randsample(vectr(:),1);
                            
                        end
                        
                        [storedpath] = Function_sphere_walk( mydist , ang , storedpath , iindex , cellradius); %simulate 'walking' patch
                        
                    end
                    
                    storedpath = storedpath + [ -startpositions(startat)*cellradius 0 0 ]; %move the entire sphere away from the center
                    
                    sphere_paths{sim,startat} = storedpath;
                    
                end
                
            end
            
        end
        
        if 0 %open flag if you want to generate figure or graphic of paths on adjacent spheres
            %%
            figure
            
            p1 = sphere_paths{7,1}; p2 = sphere_paths{7,2}; % grab two paths on spheres
            
            plot3(p1(:,1),p1(:,2),p1(:,3))
            hold on
            plot3(p2(:,1),p2(:,2),p2(:,3))
            
            [x,y,z] = sphere;
            
            x = cellradius*x; y = cellradius*y; z = cellradius*z;
            
            surf(x+cellradius,y,z,'FaceColor',[0.6 0.6 0.6] , 'EdgeColor' , 'none' , 'Facealpha' , 0.3);
            surf(x-cellradius,y,z,'FaceColor',[0.6 0.6 0.6] , 'EdgeColor' , 'none' , 'Facealpha' , 0.3);
            
            axis equal
            view([0 0]);
            axis off
            set(gcf ,'renderer' ,'painters') %do this or matlab will fucking fail to generate a vector graphic
            
        end
        
        if 1 %open flag if you want to analyze simulations for search efficiency
            
            enc_distance = 2; %encounter distance in microns
            timeofenc = nan(size(sphere_paths,1),1);
            
            for cellind = 1:size(sphere_paths,1)
                
                dists = Function_find_dist(sphere_paths{cellind,1},sphere_paths{cellind,2});
                
                if ~isempty(find(dists<enc_distance,1,'first'))
                    
                    timeofenc(cellind) = find(dists<enc_distance,1,'first');
                    
                else
                    
                    timeofenc(cellind) = Inf; %insert an arbitrarily high number
                    
                end
                
            end
            
            subplot(2,2,1)
            
            hold on
            
            dx = 0.1;
            
            if simulationtype == 1 || simulationtype == 2
            
            myang = [-pi:dx:pi];
            funct = @(x,amp,sigma,uni) amp*exp( -(x).^2/(2*sigma^2) ) + uni;
            normconst = 1/integral(@(x)funct(x,amp,parameter,uni),-pi,pi);
            %normconst = 1/( Ga( amp, sigma , uni , mya(1:(end-1))*dx ) + 0.5*dx*abs( Ga(amp,sigma,uni,mya(1:(end-1))) - Ga(amp,sigma,uni,mya(2:(end)))  ) ); %numerically integrate this function
            angdist = (funct(myang,amp,parameter,uni))*normconst;
            plot(myang,angdist)
            
            
            elseif simulationtype == 3
                
                histogram(diffangs(:) , 'DisplayStyle' , 'stairs' , 'Normalization' , 'probability');
            
            end
            
            xlim([-pi,pi])
            ylim([0 Inf])
            axis square
            
            subplot(2,2,2)
            hold on
            mycdf = cdfplot(timeofenc*timeinterval);
            xlim([0 numbsteps*timeinterval])
            xlabel('Time (min)'); ylabel('Frequency');
            axis square
            
            subplot(2,2,3)
            hold on
            timeat150 = mycdf.XData;
            timeat150(timeat150==Inf) = [];
            percentmeet = mycdf.YData(numel(timeat150));
            
            if simulationtype == 1
            
            plot(parameter,percentmeet,'.');
            xlabel('Sigma'); ylabel('Frequency');
            axis square
            
            end
        end
        
        %%
        if 0 % calculate diffusion constant of your moving particle on a sphere
            
            squaresum = zeros(size(sphere_paths{1},1),1); %store all summed squares here
            
            for cellind = 1:size(sphere_paths,1)
                
                mypath = sphere_paths{cellind,1}; mypath = mypath + [-cellradius 0 0];
                disttraveled = zeros(size(mypath,1),1); %store greatest circle distances traveled from initial location here
                
                for pointtime = 1:size(mypath,1)
                    
                    disttraveled(pointtime) = cellradius * atan2(norm(cross(mypath(1,:),mypath(pointtime,:))),dot(mypath(1,:),mypath(pointtime,:)));
                    
                    %disttraveled(pointtime) = Function_find_dist( mypath(1,:) , mypath(pointtime,:)     );
                    
                end
                
                squaresum = squaresum + disttraveled.^2;
                
            end
            
            extractedmsd = squaresum/size(sphere_paths,1);
            timeaxis = [(1:numbsteps) - 1 ]'*timeinterval;
            
            figure
            fitlimit = 8;
            
            plot(timeaxis,extractedmsd,'.')
            xlim([0 fitlimit*timeinterval]);
            
            para = polyfit(timeaxis(1:fitlimit),extractedmsd(1:fitlimit),1);
            
            disp(para(1)/4);
            
            hold on
            
            plot(timeaxis,timeaxis*para(1) + para(2))
        end
        
        
        if 0 %plotting persistence from simulated storedpath and comparing with experimental data
            
            allangles = []; %store all angles here
            diffvects  = storedpath(2:end,:) - storedpath(1:end-1,:);
            
            allrs = sqrt(sum(diffvects.^2,2));
            
            for j = 2:size(diffvects,1)
                
                diffang = atan2(norm(cross(diffvects(j,:),diffvects(j-1,:))), dot(diffvects(j,:),diffvects(j-1,:))); %find difference in angle
                allangles = [allangles, diffang];
                
            end
            subplot(1,2,1)
            hold on
            cdfplot(allangles);
            
            subplot(1,2,2)
            hold on
            cdfplot(allrs);
            
        end
    end
    
    if simulationtype==3
        legend(myparameters)
    end
    
end

%% See what search efficiency is when you sample from different probability density functions of step lengths
%                                                                      XXX
%                                                                    XXX
%                                                     +           XXX
%      XXXXXXXXXXXXX           XXXXXXXXX              |          XX                     XXXXXXX         XXXX                 XX X   XXX     XXXXX
%    XXX           XX       XXX        XXX            |         X                      XX      X       X   XX   XXX  XXXXX   X  XX XX X    X     X
%  X X              XX     XX            XXX          |        X                      XX               X    X   X X XX   X   X   XXX   X   X      X
% XX                 XX    X               X          |      XX        XXX            X                X    X   X XXX    X   X    X    X   X      X
% X             XXXX   X  XX               XXX        |     XX      XXXX             X                 X     X  X  X     X   X    X    X   X       X
% X             XXXXX  X  X      XX  XXXX    X        |    XX     XXX                X     XXXXXXX     XXXXXXX  X  X     X   X    X    X   XXXXXXXXX
% X   XXX      XXXXX  X   X    XXXXX X   X   XX       |   XX   XXX                   X        XX X     X     X  X        X   X         X   X       X
% XXXX   XXXXXX       X   XX    XXX XX   X    X       |  XX  XX         XXX          X        XX X    X      X  X       XX   X         X   X       X
%  XX                 X    X           XXX   XX       | X  XX      XXXXX             XXX    XXX  X    X      X  X       X    XX        X   X       X
%  XX               XXX    XX        XX     XX        |X X X   XXXX                    XXXXX     X    X      X
%    XX            XX       XXX      X    XX          |XXXXX
%     XXXXXXXXXXXXX           XXXXXXXXXXXX            +-----------------------+
%                                                                time


if simulate_search_efficiency_gamma_dist == 1
    
    cellradius = 2.5; %microns
    startpositions = [ 1 , -1 ]; %parameters that determine two positions where you want the patches to start
    numbsteps = 100; %number of simulations steps you want this to take
    numbsims = 5000; %number of cells you want to simulate
    timeinterval = 1.5; %minutes
    
    load('InVivo.mat');
    convertedrs = real(cellradius*acos(1 - allrs.^2/(2*cellradius^2))); % since these are all straight lines, convert to how far a point would have traveled on a sphere of some radius
    gammamean = mean(convertedrs(:)); %mean you will use for gamma distribution
    
    allpercentmeet = []; allsuccessenc = []; %store enconter data here
    gshapenum = 0;
    
    for gshape = [ -1 -2 0 0.3 1.1 2 4  ] %when gshape == 0, use experimentally measured values
        
        %for gshape = [ -2 ] %when gshape == 0, use experimentally measure values
        gshapenum  = gshapenum  + 1;
        
        sphere_paths = cell(numbsims,numel(startpositions));
        
        %convertedrs = allrs;
        if 0 %open flag if you want to simulate stuff
            alldist =  NaN(numbsims*numbsteps*2,1); distii = 0; %this is so you can generate histograms of distributions you sample from
            for sim = 1:numbsims
                
                for startat = 1:numel(startpositions) %start with patch at either extreme of the x axis
                    
                    %generate a point on the back half of spherical cell
                    
                    ele = pi*rand() + startpositions(startat)*pi/2;
                    az = asin(2*rand()-1);
                    [Xc,Yc,Zc] = sph2cart(ele,az,cellradius);
                    
                    storedpath(1,:) = [ Xc,Yc,Zc ];
                    
                    for iindex = 2:numbsteps
                        
                        amp = 1;
                        sigma = 0.9;
                        uni = 0.2;
                        
                        if gshape == 0
                            
                            %simulate experimental cell's movement with both r and theta
                            ang = Function_sample_normal_uniform( amp , sigma , uni ); % amp = 1, sigma = 0.9, and uni = 0.2 seems to match turning angle distribution seen in experiments
                            mydist = randsample(convertedrs,1);
                            %mydist = randsample(convertedrs(convertedrs<1.5),1); %this is to see how searching works when you truncate higher step lengths
                            
                        elseif gshape == -1
                            
                            %simulate movement with r being just the mean
                            
                            ang = Function_sample_normal_uniform( amp , sigma , uni ); % amp = 1, sigma = 0.9, and uni = 0.2 seems to match turning angle distribution seen in experiments
                            mydist = gammamean;
                            
                        elseif gshape == -2
                            
                            ang = Function_sample_normal_uniform( amp , sigma , uni ); % amp = 1, sigma = 0.9, and uni = 0.2 seems to match turning angle distribution seen in experiments
                            mydist = randsample(convertedrs(convertedrs<1.5),1); %this is to see how searching works when you truncate higher step lengths
                            
                        else
                            
                            gscale = gammamean/gshape;
                            ang = Function_sample_normal_uniform( amp , sigma , uni );
                            mydist = gamrnd(gshape,gscale);
                            
                            distii = distii+1 ;alldist(distii) = mydist;
                            
                        end
                        
                        [storedpath] = Function_sphere_walk( mydist , ang , storedpath , iindex , cellradius); %simulate 'walking' patch
                        
                    end
                    
                    
                    storedpath = storedpath + [ -startpositions(startat)*cellradius 0 0 ]; %move the entire sphere away from the center
                    
                    sphere_paths{sim,startat} = storedpath;
                    
                end
                
            end
            
            alldist(isnan(alldist))=[];
            
            save( [ 'Gammastuff_', num2str(gshape) , '.mat' ] , 'alldist' , 'sphere_paths' ,'gshape')
        end
        
        %%
        if 1 %open flag if you want to analyze simulations for search efficiency
            
            load( [ 'Gammastuff_', num2str(gshape) , '.mat' ] )
            
            enc_distance = 2; %encounter distance in microns
            timeofenc = nan(size(sphere_paths,1),1);
            residencetime = nan(size(sphere_paths,1),1);
            allresidencetimes = [];
            
            for cellind = 1:size(sphere_paths,1)
                
                dists = Function_find_dist(sphere_paths{cellind,1},sphere_paths{cellind,2});
                
                if ~isempty(find(dists<enc_distance,1,'first'))
                    
                    timeofenc(cellind) = find(dists<enc_distance,1,'first'); %find out when they first come within encounter distance
                    
                else
                    
                    timeofenc(cellind) = Inf; %insert an arbitrarily high number
                    
                end
                
                %resttime = sum(dists<enc_distance);
                %if resttime == 0, resttime = NaN; end
                %residencetime(cellind) = 1.5*resttime; %amount of time patches spend near each other
                
                resttime = dists<enc_distance;
                begt = strfind(resttime',[0 1]);
                endt = strfind(resttime',[1 0]);
                if numel(begt)>numel(endt),begt(end) =[]; end
                residencetime(cellind) = 1.5*mean(endt-begt);
                
                allresidencetimes = [ allresidencetimes , 1.5*(endt-begt) ];
                
            end
            
            subplot(4,7,gshapenum) %cdfs of when patches meet up
            hold on
            mycdf = cdfplot(timeofenc*timeinterval);
            xlim([0 numbsteps*timeinterval])
            xlabel('Time (min)'); ylabel('CF');
            axis square
            
            subplot(4, 7 ,[ 8 : 14 ] ) % bargraphs of histogram
            
            hold on
            grid on
            timeat150 = mycdf.XData;
            timeat150(timeat150==Inf) = [];
            percentmeet = mycdf.YData(numel(timeat150));
            
            bar(gshapenum,percentmeet)
            
            xlim([0.6 ,7.4])
            
            subplot(4,7,gshapenum+14) %plot histograms of how distributions look
            
            dx = 0.1;
            
            hold on
            
            if gshape == -1
                
                plot([gammamean,gammamean],[0 4])
                
            elseif gshape == 0
                
                [ histvalues,histedges ] = histcounts(convertedrs,'Normalization','pdf','Binwidth',dx);
                plot( histedges(2:end)-dx/2 , histvalues) %the -dx/2 is to make it align with theoretical distribution
                
            elseif gshape == -2
                
                [ histvalues,histedges ] = histcounts(convertedrs(convertedrs<1.5),'Normalization','pdf','Binwidth',dx);
                plot( histedges(2:end)-dx/2 , histvalues) %the -dx/2 is to make it align with theoretical distribution
                
            else
                
                [ histvalues,histedges ] = histcounts(alldist,'Normalization','pdf','Binwidth',dx);
                plot( histedges(2:end)-dx/2 , histvalues) %the -dx/2 is to make it align with theoretical distribution
                
            end
            
            xlim([0 2])
            xlabel('Step length')
            ylabel('Frequency')
            axis square
            ylim([0 5])
            
            subplot(4,7,22)
            
            hold on
            
            if gshape == -1
                
                plot( 0 , percentmeet , '.' , 'MarkerSize' , 20  )
                
            elseif gshape == 0
                
                plot( std(convertedrs) , percentmeet , '.' , 'MarkerSize' , 20  )
                
            elseif gshape == -2
                
                plot( std(convertedrs(convertedrs<1.5)) , percentmeet , '.' , 'MarkerSize' , 20  )
                
            else
                
                plot( std(alldist) , percentmeet , '.' , 'MarkerSize' , 20  )
                
            end
            
            axis square
            xlabel('Std'); ylabel('% meet'); ylim([0 Inf])
            
            %{
            
            subplot(2,3,1)
            hold on
            dx=1.5;
            [ histvalues,histedges ] = histcounts(allresidencetimes,'Normalization','pdf','Binwidth',dx);
                plot( histedges(2:end)-dx/2 , histvalues) %the -dx/2 is to make it align with theoretical distribution
            
            subplot(2,3,2) %cdf plot for time at which patches get within a certain distance of each other
            hold on
            mycdf = cdfplot(timeofenc*timeinterval);
            xlim([0 numbsteps*timeinterval])
            xlabel('Time (min)'); ylabel('Frequency');
            axis square
            
            subplot(2,3,3)
            hold on
            timeat150 = mycdf.XData;
            timeat150(timeat150==Inf) = [];
            percentmeet = mycdf.YData(numel(timeat150));
            plot(gshape,percentmeet,'.','MarkerSize',20);
            xlabel('Sigma'); ylabel('Frequency');
            axis square
            
            subplot(2,3,4) %plot histogram of different distributions
            dx = 0.1;
            
            hold on
            if gshape == -1
                
                plot([gammamean,gammamean],[0 4])
                
            elseif gshape == 0
                
                [ histvalues,histedges ] = histcounts(convertedrs,'Normalization','pdf','Binwidth',dx);
                plot( histedges(2:end)-dx/2 , histvalues) %the -dx/2 is to make it align with theoretical distribution
                
            else
                
                [ histvalues,histedges ] = histcounts(alldist,'Normalization','pdf','Binwidth',dx);
                plot( histedges(2:end)-dx/2 , histvalues) %the -dx/2 is to make it align with theoretical distribution
                
            end
           
            
            xlim([0 2.5])
            axis square
            
            hold on
            
            subplot(2,3,5)
            
            hold on
            yyaxis left
            plot(gshape,percentmeet,'.','MarkerSize',20);
            yyaxis right
            
            plot(gshape, mean(allresidencetimes) ,'.','MarkerSize',20)
            
            allpercentmeet = [ allpercentmeet , percentmeet ];
            allsuccessenc = [ allsuccessenc , mean(allresidencetimes) ];
            axis square
            
            %}
            
        end
        
    end
    
    
end

%% Generate illustration of cells searching for each other by making their patches wander

if create_figure_of_cells_searching == 1
    
    cellradius = 2.5;
    
    load('Gammastuff_0.mat')
    cellpair = 9;
    p1 = sphere_paths{cellpair,1}; p2 = sphere_paths{cellpair,2}; % grab two paths on spheres
    
    plot3(p1(:,1),p1(:,2),p1(:,3))
    hold on
    plot3(p2(:,1),p2(:,2),p2(:,3))
    
    [x,y,z] = sphere;
    
    x = cellradius*x; y = cellradius*y; z = cellradius*z;
    
    surf(x+cellradius,y,z,'FaceColor',[0.6 0.6 0.6] , 'EdgeColor' , 'none' , 'Facealpha' , 0.3);
    surf(x-cellradius,y,z,'FaceColor',[0.6 0.6 0.6] , 'EdgeColor' , 'none' , 'Facealpha' , 0.3);
    
    axis equal
    view([0 0]);
    axis off
    set(gcf ,'renderer' ,'painters') %do this or matlab will fucking fail to generate a vector graphic
    
    
    %myp = sphere_paths{cellpair,2};
    
    %plotOnSphere(myp(:,1) , myp(:,2) , myp(:,3) )
    
    
    
end


%% see how much patch characteristics change over the course of a simulation
%
%                      Observe how patch changes
%
%           XXXXX
%           X   X
%          XX   X                    XX                      X
%          X    X                     XXX                   X  X
%          X    XX          XXXXXXXXXXXXXXX         XXX    XX   X
%         XX     X                    XXX          XX X    X    X
%        XX      XX                  XX           XX  XX   X    X
%      XXX        XX                              X    XXXXX    XX
%    XXX           XXXX                          X               XX
%    X                XX                        X                 XX
% X XXXXXXXXXXXXXXXXXXXXXX                    XXXXXXXXXXXXXXXXXXXXXXXXX



if see_patch_variability == 1
    
    foldername ='C:\Users\dg144\Desktop\Projects\Simulations\20190113_SETMAWP_Treadmill\OUT\' ;
    
    
    for filename = {'20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024','20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_GAP_090','20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_LocalGAP_thalf_1',...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_LocalGAP_thalf_10','20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024_LocalGAP_thalf_30_GAP_2'}
        
        
        max42s = nan(400*90,1); %peak max stored here
        sum42s = nan(400*90,1); %total 42 on membrane stored here
        enternum = 0; %index for entries
        
        for ii = 1:400
            
            openfile = fullfile(foldername,char(filename),num2str(ii),'AllProteins_Pher_10nM.mat');
            
            clear Total_42T_store %clear this from memory
            
            if exist(openfile,'file')
                
                load(openfile,'Total_42T_store')
                
                for jj = 1:size(Total_42T_store,3)
                    
                    enternum = enternum+1;
                    max42s(enternum) = max(max(Total_42T_store(:,:,jj)));
                    sum42s(enternum) = sum(sum(Total_42T_store(:,:,jj)));
                    
                end
                
                
            end
            
        end
        
        
        subplot(2,2,1)
        hold on
        histdata = histogram(max42s,'BinWidth',3,'Normalization','probability','DisplayStyle','stairs');
        xlabel('GTP-Cdc42 \mu M')
        axis square
        
        set(gcf ,'renderer' ,'painters') %do this or matlab will fucking fail to generate a vector graphic
        
        subplot(2,2,2)
        hold on
        plot(histdata.BinEdges(2:end),histdata.Values(1:end))
        xlabel('GTP-Cdc42 \mu M')
        axis square
        
        subplot(2,2,3)
        hold on
        histdata = histogram(sum42s,'BinWidth',200,'Normalization','probability','DisplayStyle','stairs');
        xlabel('GTP-Cdc42 \mu M')
        axis square
        
        set(gcf ,'renderer' ,'painters') %do this or matlab will fucking fail to generate a vector graphic
        
        subplot(2,2,4)
        hold on
        plot(histdata.BinEdges(2:end),histdata.Values(1:end))
        xlabel('GTP-Cdc42 \mu M')
        axis square
        
        
        
    end
    
    legend('Gaussian=1','Global GAP = 0.090','Local GAP thalf = 1s amp = 1/s','Local GAP thalf = 10s amp = 1/s','Local GAP thalf = 30s amp = 2/s')
    
end

%% model point spread function

%modeling point spread of microscope

if model_microscope_point_spread == 1
    
    collectdata = 0; %decided if you want to first collect data for your PSFs
    fitto2dgauss = 0;
    
    if collectdata == 1
        
        myROIs = {};
        cellnumber = size(myROIs,2);
        radius = 6; %radius of circle around point you choose
        filename = 'Z:\Users\Debraj\Movies and Data\2017\20190419_PSF\MAX_PSF_100x_Spinning_Disk_7.TIF'; %green channel
        im = imread(filename);
        
        while(1)
            imagesc(im); axis square; %show image
            cellnumber = cellnumber + 1; %increment the cell index
            donewithcells = input('Done marking all cells?','s'); % when you are done with marking all cells, terminate the loop
            
            if donewithcells == 'y'
                break
            end
            
            if ~isempty(myROIs) %this is just to draw circles around all points that you have already marked
                for filenametoload = 1:size(myROIs,2)
                    drawcircle('Position',myROIs{filenametoload},'Radius',radius);
                end
            end
            
            circ = drawpoint;%circ = drawcircle('Center',[100,100],'Radius',10);
            
            while(1) %pause code to make changes to your ROI, give user option to proceed once they are satisfied
                m = input('Happy with ROI? ','s');
                if m == 'y'
                    break
                end
            end
            
            %save your ROI here, maybe within a cell matrix?
            
            myROIs{cellnumber} = circ.Position;
            
        end
        
        save('PSF_ROIs','myROIs');
        
    end
    
    if  fitto2dgauss == 1
        
        load('PSF_ROIs.mat');
        filename = 'Z:\Users\Debraj\Movies and Data\2017\20190419_PSF\MAX_PSF_100x_Spinning_Disk_7.TIF'; %green channel
        im = imread(filename);
        allsigmas = [];
        for filenametoload = 1:size(myROIs,2)
            
            centx = round(myROIs{filenametoload}(1)); centy = round(myROIs{filenametoload}(2));
            
            smallsquare = im( centy-10 : centy+10 , centx-10 : centx+10  );
            smallsquare = double(smallsquare);
            
            [myx , myy ] = meshgrid( 1:size(smallsquare,1) , 1:size(smallsquare,1) );
            
            [fitresult, gof] = createFit(myx, myy, smallsquare,filenametoload);
            
            allsigmas = [ allsigmas , fitresult.sig ];
            PSFsigma = mean(allsigmas(:));
            
        end
        
    end
    
    %% Now estimate what your simulation polarisome looks like and apply convolution to your simulation images.
    
    % reference - http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/VELDHUIZEN/node10.html
    PSFsigma = 1.0449; %derived from analysis I did above
    [myx , myy ] = meshgrid( 1:21 , 1:21) ;
    PSF = exp( -(myx-11).^2/(2*PSFsigma ^2) -(myy-11).^2/(2*PSFsigma ^2) ); %generate PSF
    simpixtomic = 0.0886; exppixtomic = 0.1386; %pixel to micron conversion for simulation and experiment
    
    PSF = imresize(PSF, exppixtomic/simpixtomic );
    
    %% plot how actin cable clustering looks in your simulations
    
    foldername = 'C:\Users\dg144\Desktop\Projects\Simulations\20190113_SETMAWP_Treadmill\OUT\';
    for filenametoload = {'20190113_SETMAWP_Treadmill_0_0',...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_1_actinvar_024',...
            '20190113_G_SETMAWP_Treadmill_0_0_k_gauss_5_actinvar_024'}
        
        myimage = zeros(100,100);
        
        for ii = 1:50
            
            load([foldername,char(filenametoload),'\',num2str(ii),'\AllProteins_Pher_10nM.mat'],'storecables','c42cen');
            
            for jj = 1:540
                
                %put everything near the center
                cable = storecables(jj,:); cable(cable<=0) = []; %remove empty cable positions on the membrane (since -1 means no cables)
                recenter_by = [50.5,50.5] - c42cen(jj,:); % coarse recentering
                recenter_by = repmat(recenter_by , [1 , numel(cable)/2]);
                cable = mod(cable + recenter_by,100+10^-12) ; %move cables to center of patch. mod(i-1,100)+1 to get values between 1 and 100.
                cablecentroid = [ mean(cable(1:2:end)) , mean(cable(2:2:end))] ; %now find actual center of cables
                recenter_by = [50.5,50.5] - cablecentroid; %fine recentering
                recenter_by = repmat(recenter_by , [1 , numel(cable)/2]);
                cable = cable + recenter_by;
                
                cableimage = fillcables(zeros(100,100),cable);
                myimage = myimage + cableimage; % build image from cable locations
                
                %pause(0.1)
                %imagesc(cableimage)
                %drawnow
                
            end
            
        end
        
        load([foldername,'Avg_RecA_Cdc42T_Profile_20190113_SETMAWP_Treadmill_0_20.mat'],'avg42');
        
        figure
        subplot(2,2,1)
        
        imagesc(myimage)
        title('Polarisome')
        colormap(flipud(gray))
        axis square
        %axis([38 62 38 62])
        
        
        
        subplot(2,2,2)
        % convolution step
        
        convimpol = conv2(PSF,myimage); convimpol = convimpol(17:116,17:116);
        convimpol = convimpol( 38:62 , 38:62 ); %extract just the patch
        convimpol = convimpol/(sum(convimpol(:)));
        
        imagesc(convimpol);
        title('Blurred polarisome')
        axis square
        colormap(flipud(gray))
        %axis([38 62 38 62])
        colorbar
        caxis([0 10e-3])
        
        subplot(2,2,3)
        
        imagesc(avg42);
        title('Patch')
        axis square
        colormap(flipud(gray))
        %axis([38 62 38 62])
        subplot(2,2,4)
        
        convim42 = conv2(PSF,avg42); convim42 = convim42(17:116,17:116);
        convim42 = convim42( 38:62 , 38:62 ); %extract just the patch
        convim42 = convim42/(sum(convim42(:)));
        
        imagesc(convim42);
        title('Blurred patch')
        axis square
        colormap(flipud(gray))
        %axis([38 62 38 62])
        colorbar
        caxis([0 10e-3])
        
        
        figure
        
        plot([1:round(size(convimpol,1))]*0.0886,convimpol(round(size(convimpol,1)/2) ,:))
        hold on
        plot([1:round(size(convimpol,1))]*0.0886,convim42(round(size(convimpol,1)/2),:))
        xlim([0 round(size(convimpol,1))*0.0886])
        axis square
        ylabel('Intensity (a.u.)')
        ylim([ 0 , 10e-3 ])
        %
    end
    
    
    
    %% testing
    
    %{
    myimage = zeros(100,100);
    
    for i = 1:10, myimage( randsample([45:55],1) , randsample([45:55],1)) = 1; end
    
    convim = conv2(PSF,myimage); convim = convim(17:116,17:116);
    deconvim = deconvwnr(convim,PSF,0);
    
    figure
    
    subplot(1,3,1)
    imagesc(myimage);
    title('Points')
    axis square
    colormap gray
    
    subplot(1,3,2)
    imagesc(convim);
    title('Crappy image')
    axis square
    colormap gray
    
   subplot(1,3,3)
   imagesc(deconvim)
   title('Deconvolved')
   axis square
   colormap gray
    %}
    
end

if sample_profiles_of_polarisome_and_patch == 1
    
    %% maybe put in some kind of condition so you stop the script once you are done with marking ROIs on cells
    
    pixtomic = 0.1386;
    probenames = 'Spa2_Bem1';
    
    filename1 = fullfile('Dataofpolarisomeandpatch',[probenames,'.tif']); %red channel
    im1 = imread(filename1,1);
    filename2 = fullfile('Dataofpolarisomeandpatch',[probenames,'.tif']); %green channel
    im2 = imread(filename2,2);
    
    im1adj = imadjust(im1); %adjust stuff so it is easier to see
    im2adj = imadjust(im2);
    dispimg = cat(3,im1adj,im2adj,zeros(size(im1,1)));
    
    for cc = 1
        
        image(dispimg); %show image
        
        %now draw your polygon
        
        cellbound = drawpolyline;%circ = drawcircle('Center',[100,100],'Radius',10);
        
        while(1) %pause code to make changes to your ROI, give user option to proceed once they are satisfied
            m = input('Happy with ROI? ','s');
            if m == 'y'
                break
            end
        end
        
        linedrawn = cellbound.Position;
        
        %extract profile for red channel
        
        numpix = 50; %how much sampling you want there to be in the line
        
        myimage = im1;
        [distance,avgprofile] = avgline(linedrawn,myimage,numpix);
        plot(distance*pixtomic,(avgprofile-min(avgprofile))/(max(avgprofile)-min(avgprofile)),'r');
        
        %extract profile for green channel
        
        myimage = im2;
        [distance,avgprofile] = avgline(linedrawn,myimage,numpix);
        hold on
        plot(distance*pixtomic,(avgprofile-min(avgprofile))/(max(avgprofile)-min(avgprofile)),'g');
        
        xlabel('Distance $\mu$m','Interpreter','Latex')
        ylabel('Fluoresence (a.u.)')
        xlim([min(distance*pixtomic),max(distance*pixtomic)])
        axis square
        
        save(fullfile( pwd , 'Dataofpolarisomeandpatch' , ['Boundary_',probenames]),'linedrawn');
        saveas(gcf,fullfile( pwd , 'Dataofpolarisomeandpatch' , ['Profiles_',probenames,'.pdf']));
        
        
    end
    
end

if make_figure_bursty_deliv == 1
    
   
    
end



%% Functions

%% function to generate curved lines on spherical surface




function plotOnSphere(x,y,z)

%// Vectors representing each point
xyz = [x(:), y(:), z(:)].';  %'

%// One vector of the "first" points and one of the "next" points
v1 = xyz(:, 1:end-1);
v2 = xyz(:, 2:end);

%// Cross product between the vectors of one point and the next
cv1v2 = cross(v1, v2);

%// Compute unit vector in the plane defined by v1 and v2
v3 = normc(cross(cv1v2, v1));

%// Figure out the range of the inner angle between v1 and v2
nc = sqrt(sum(cv1v2.^2, 1));
t = atan2(nc, dot(v1, v2, 1));

%// Number of points to sample between any two points on the sphere
nPoints = 100;

%// Compute the interpolant
V = zeros([nPoints, fliplr(size(v1))]);
for k = 1:numel(t)
    T = linspace(0, t(k), 100);
    V(:,k,:) = (v1(:,k) * cos(T) + v3(:,k) * sin(T)).';    %'
end

%// Break the result out into x,y,z parts
xx = V(:,:,1);
yy = V(:,:,2);
zz = V(:,:,3);

%// Plot the lines
h = plot3(xx(:), yy(:), zz(:));
hold on

%// Plot the original data points
plot3(x,y,z, 'o', ...
    'Color', get(h, 'Color'), ...
    'Parent', get(h, 'Parent'));
end


%% function to generate averaged profile of a line on an image

function [distance,avgprofile] = avgline(linedrawn,myimage,num)

vector = [linedrawn(2:end,:)-linedrawn(1:end-1,:)];
th = cart2pol( vector(:,1) ,vector(:,2)  );

[deltaxn,deltayn] = pol2cart(th-pi/2,1);
[deltaxp,deltayp] = pol2cart(th+pi/2,1);

deltaxn = [deltaxn;deltaxn(end)]; deltayn = [deltayn;deltayn(end)];
deltaxp = [deltaxp;deltaxp(end)]; deltayp = [deltayp;deltayp(end)];

linen = linedrawn + [deltaxn,deltayn]; %lines on either side of my current line
linep = linedrawn + [deltaxp,deltayp];

[cx,cy,intensity] = improfile(myimage,linedrawn(:,1),linedrawn(:,2),num);
[cx1,cy1,intensity1] = improfile(myimage,linedrawn(:,1),linedrawn(:,2),num);
[cx2,cy2,intensity2] = improfile(myimage,linedrawn(:,1),linedrawn(:,2),num);

avgprofile = [intensity , intensity1 , intensity2 ];
avgprofile = mean(intensity,2);

% while you're at it calculate distance covered by the line you drew

vect = [cx(2:end)-cx(1:end-1) , cy(2:end)-cy(1:end-1)];

[~,r] = cart2pol(vect(:,1),vect(:,2));

distance = cumsum(r);
distance = [0; distance];

end

%% fill matrix with cables

function [myimage] = fillcables(myimage,cable)

cable = round(cable); cable(cable==0)=100; %adjust cable values

xs = cable(1:2:end); ys = cable(2:2:end);

for ii = 1:numel(xs), myimage(xs(ii),ys(ii)) = 1; end

end


%% 2d Gaussian fit to data

function [fitresult, gof] = createFit(myx, myy, smallsquare,ii)

[xData, yData, zData] = prepareSurfaceData( myx, myy, smallsquare );

% Set up fittype and options.
ft = fittype( 'a*exp( -(myx-ux)^2/(2*sig^2) -(myy-uy)^2/(2*sig^2) ) + unif', 'independent', {'myx', 'myy'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [ max(smallsquare(:)) , 1 , 100 , 9 , 9 ];
opts.Robust = 'LAR';
opts.Lower = [max(smallsquare(:))*0.7 0.8 0 0 0];
%opts.Upper = [Inf 3 Inf Inf Inf]; %set fit bounds

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
subplot(4,4,ii)
plot( fitresult, [xData, yData], zData );
% Label axes
xlabel myx
ylabel myy
zlabel smallsquare
grid on
view( -34.5, 4.1 );

end


%% simple displacement finder

function [c] = Function_find_dist(a,b)

c = sqrt(sum((a - b).^2 , 2 ));

end

%% sample from normal + uniform distribution

function [b] = Function_sample_normal_uniform( amp , sigma , uni )

upperlimit = amp*exp( -(0)^2/(2*sigma^2) ) + uni;

a = 1;
upperbound = 0;

while a > upperbound %keep sampling b till Function(b) < a
    
    
    a = rand()*upperlimit; %sample for a value between 0 and max of distribution
    b = rand()*2*pi - pi; %sample to pick point along x axis between pi and -pi
    upperbound = amp*exp( -(b).^2/(2*sigma^2) ) + uni;
    
end

end

%% Function to diffuse particle on surface of a sphere

function [storedpath ] = Function_sphere_walk(dist , ang , storedpath, stepnum, cellradius )
%Function_sphere_walk Makes a point do a random walk on a sphere of radius r. The point's step size per time step (dt) is sampled from a normal distribution.
%   storedpath contains cartesian coordinates of path. vectr is the
%   matrix or real number that determines how far a particle move in a given
%   timestep. diffangs is a matrix or real number that determines what
%   angle the spot will take wrt to its previous movement. stepnumb is what
%   step in the simulation this is in

%% first move last to points to the north pole so calculations are easy

northpole = [0 0 cellradius];
prevpoint = storedpath(stepnum-1,:);

normvect2pole = cross(prevpoint,northpole); normvect2pole = normvect2pole/norm(normvect2pole);%find normal vector to north pole

if stepnum>2, pointstomove = storedpath( (stepnum-2) : (stepnum-1) , : ); %these are the points you want to move to pole to generate where the patch will be in the next time step
else, pointstomove = storedpath(stepnum-1,:);
end

v = pointstomove; k = normvect2pole.*ones(size(pointstomove,1),3); % k is the vector you are carrying out rodriquez rotation around

[~,rottopole,~] = cart2sph(pointstomove(end,1),pointstomove(end,2),pointstomove(end,3));
rottopole = pi/2-rottopole; %so you rotate toward the north pole

pointstomove = v*cos(rottopole) + cross(k,v)*sin(rottopole) + k.*dot(k,v,2)*(1-cos(rottopole)); %apply Rodriguez transform to move point(s) to the pole

moveby = dist/cellradius; %express how much you want to move in radians

if size(pointstomove,1) == 2, [prevangle , ~, ~] = cart2sph(pointstomove(1,1),pointstomove(1,2),pointstomove(1,3)); %you can only have prevangle wrt to previous point
else prevangle = rand()*2*pi; %if there is no previous point, simply pick an angle at random
end

%% now make point 'walk'

newpoint = [prevangle+pi,pi/2-moveby,cellradius];
[newx , newy, newz] = sph2cart( newpoint(1) , newpoint(2) , newpoint(3) );
v = [newx , newy, newz]; k = [ 0 0 1 ];
rotaboutpole = ang; %angle you want it to move wrt to previous angle

newpoint = v*cos(rotaboutpole) + cross(k,v)*sin(rotaboutpole) + k.*dot(k,v,2)*(1-cos(rotaboutpole)); %apply Rodriguez transform to rotate point about the pole

%% now move the 'walked' point back to where it belongs

rottopole = -rottopole;
v = newpoint; k = normvect2pole;
newpoint = v*cos(rottopole) + cross(k,v)*sin(rottopole) + k.*dot(k,v,2)*(1-cos(rottopole)); %rotate newpoint back to where it should be

storedpath(stepnum,:) = newpoint;

end


