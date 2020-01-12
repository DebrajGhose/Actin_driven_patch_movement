set(0,'DefaultTextFontSize',18);
set(0,'DefaultAxesFontSize',18);
[status currdir] = system('pwd'); %to use as an automatic label for some graphs that are going to be generated for many diff. setups
N0=0; %starting number for the directory name
Nruns=20; %number of the runs\trajectories to read
dti=1;      % integration step (s)
dtstring='1';
dt=30/60; %delta t in min from the RUN_.....
N=100;
ddx=5/N*sqrt(pi);

xrs=[];yrs=[];

counti=0;
%% color selection
orange=[255 150 0]./255;
blue = [21 161 255]./255;
chartreuse2=[118 238 0]./255;
dodgerBlue4=[16 78 139]./255;
aquamarine3=[102 205 170]./255;
RoyalBlue=[65 105 225]./255;
DeepPink2 = [238 18 137]./255;
MediumPurple=[147 112 219]./255;

[X,Y]=meshgrid(1:N,1:N);
rdist=sqrt((X-N/2).^2+(Y-N/2).^2);
xrs={};yrs={};Nfs=[];
mkdir('./plots');

%% Calculation of the patch geometric center evolution in time; to be stored in mcoord = array(number of runs,number of frames, 2); 2=>[xind, yind] indexes of the patch geom center.
for ii=1:Nruns %in case we are averaging over many runs. 
    for j=1:5
        fname=['./' num2str(ii+N0) '/cdc42t_reca_weight_centers_' num2str(j) '.mat'];   
        if( exist(fname,'file')==2)
            load(fname);
            xri=[];yri=[];
            xri=c42cen(:,1);
            yri=c42cen(:,2);
            Nfs=[Nfs; size(c42cen,1)]; 
    %             figure(34);clf; hold on;
    %             scatter(xri,yri,70,dodgerBlue4,'filled');
    %             scatter(xri(1),yri(1),'green','fill');
    %             scatter(xri(end),yri(end),'red','fill');
    %             xlim([1 100]);ylim([1 100]);
    %             box on;
    %             fnamei=['./graphs/trace_dir' num2str(ii) '-run' num2str(j) '-' num2str(counti) '.png'];
    %             saveas(figure(34),fnamei,'png');
            xrs{end+1}=xri;
            yrs{end+1}=yri;
        end        
    end
end
Nruns=length(Nfs);
display(['---> Number of trajectories to analyze: ' num2str(Nruns)]);


%% M.S.D. calculation: take into account periodic boundary conditions(pbc).
msdx=[];msdy=[];
xave=[];yave=[];x2m=[];y2m=[];msd2x=[];msd2y=[];

%Nstops=[1, 5, 7, 10, 15, 20, 30, 40, 50 ,60]./dt;
Nstops = [1 (10:10:180)./dt];
for i=1:length(Nstops)
    xave_i=[];yave_i=[];
    for k=1:Nruns
        if (Nfs(k)>=Nstops(i))
            dx=xrs{k}(Nstops(i))-xrs{k}(1);
%             if(abs(dx)>N/2)           
%                 images=[dx (dx+N) (dx-N)];
%                 [min_i,ind_i]=min(abs(images));
%                 dx=images(ind_i);             
%             end
            dy=yrs{k}(Nstops(i))-yrs{k}(1);
%             if(abs(dy)>N/2)               
%                 images=[dy (dy+N) (dy-N)];
%                 [min_i,ind_i]=min(abs(images));
%                 dx=images(ind_i);                
%             end
%             if(dri>N/2)
%                 minri=10000000;
%                 for dispi=pbc %possible periodic images for x and y
%                     for dispj=pbc
%                         ralt=sqrt((dx+dispi).^2+(dy+dispj).^2);
%                         if(ralt<minri) 
%                             minri=ralt;
%                             dx=dx+dispi;
%                             dy=dy+dispj;
%                         end
%                     end
%                 end
%                 dri=minri;
%             end
            xave_i=[xave_i dx*ddx];
            yave_i=[yave_i dy*ddx];       
        end
    end    
    msdx=[msdx mean(xave_i.^2)];
    msdy=[msdy mean(yave_i.^2)];
    msd2x=[msd2x std(xave_i.^2)./sqrt(length(xave_i))];
    msd2y=[msd2y std(yave_i.^2)./sqrt(length(xave_i))];

    xave=[xave mean(xave_i)];
    yave=[yave mean(yave_i)];
    x2m=[x2m std(xave_i)./sqrt(length(xave_i))];
    y2m=[y2m std(yave_i)./sqrt(length(yave_i))];
end

skip=1;
f210=figure(210);
clf; hold on;
h1=plot(Nstops(1:skip:end).*dt, xave(1:skip:end),'o-', 'Color',RoyalBlue,'LineWidth',6);
errorbar(Nstops(1:skip:end).*dt, xave(1:skip:end), x2m(1:skip:end),'Color',RoyalBlue,'LineWidth',2);
h2=plot(Nstops(1:skip:end).*dt, yave(1:skip:end),'o-', 'Color',aquamarine3,'LineWidth',6);
errorbar(Nstops(1:skip:end).*dt, yave(1:skip:end),y2m(1:skip:end), 'Color',aquamarine3,'LineWidth',2);
%plot(Nstops.*dt, rave,'o-', 'Color',DeepPink2,'LineWidth',6);
xlabel('Time (min)');
ylabel('<x(t) - x(t-s)> (\mum)');
h=legend([h1 h2],'<\Delta x>','<\Delta y>');%,'<\Delta r>');
set(h,'Location','NorthWest');
title([' Number of analyzed trajectories: ' num2str(Nruns)]);
xlim([0 max(Nstops)*dt]);
xlim([0 60]);
ylim([-0.1 2.5]);grid on;
fname=['./plots/dx_dy_runs' num2str(Nruns) '_1-.png'];
saveas(f210,fname,'png');
xlim([0 max(Nstops)*dt]);
ylim([-0.1 5]);grid on;
set(gca,'YTick',0:0.5:5,'FontSize',24);
fname=['./plots/dx_dy_runs' num2str(Nruns) '_2-.png'];
saveas(f210,fname,'png');

f211=figure(211);
clf; hold on;
h1=plot(Nstops(1:skip:end).*dt, msdx(1:skip:end),'o-', 'Color',RoyalBlue,'LineWidth',6);
errorbar(Nstops(1:skip:end).*dt, msdx(1:skip:end), msd2x(1:skip:end),'Color',RoyalBlue,'LineWidth',2);
h2=plot(Nstops(1:skip:end).*dt, msdy(1:skip:end),'o-', 'Color',aquamarine3,'LineWidth',6);
errorbar(Nstops(1:skip:end).*dt, msdy(1:skip:end), msd2y(1:skip:end),'Color',aquamarine3,'LineWidth',2);
%plot(Nstops.*dt, rave,'o-', 'Color',DeepPink2,'LineWidth',6);
xlabel('Time (min)');
ylabel('<(x(t) - x(t-s))^2> (\mum^2)');
h=legend([h1 h2],'<\Delta x^2>','<\Delta y^2>');%,'<\Delta r>');
set(h,'Location','NorthWest');
ylim([0 9]);xlim([1 180]);
grid on;
fname=['./plots/sqdisp-xy-0nM_equil_' num2str(Nruns) 'runs_2-.png'];
saveas(f211,fname,'png');
   
save('deltaxy_frominitpos_1.mat', 'xave','yave','x2m','y2m','msd2x','msd2y','msd','rave','mstd', 'Nstops'); 