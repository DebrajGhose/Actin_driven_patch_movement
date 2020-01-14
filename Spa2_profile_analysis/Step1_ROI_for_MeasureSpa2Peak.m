%% figure out what the size of a patch is
clear all
close all

%% maybe put in some kind of condition so you stop the script once you are done with marking ROIs on cells

myROIs = {};
cellnumber = size(myROIs,2);
radius = 10; %radius of circle around point you choose 
filename1 = 'C1-20160206_PDGGTexp3_dly18172_dly18274_alpha0nM_488_561_merged.tif'; %green channel
im1 = imread(filename1);
filename2 = 'C2-20160206_PDGGTexp3_dly18172_dly18274_alpha0nM_488_561_merged.tif'; %red channel
im2 = imread(filename2);

im1 = imadjust(im1); %adjust stuff so it is easier to see
im2 = imadjust(im2);
dispimg = cat(3,im2,im1,zeros(size(im1,1)));

while(1)
    image(dispimg); %show image 
    cellnumber = cellnumber + 1; %increment the cell index
    donewithcells = input('Done marking all cells?','s'); % when you are done with marking all cells, terminate the loop
    
    if donewithcells == 'y'
        break
    end
    
    if ~isempty(myROIs) %this is just to draw circles around all points that you have already marked
    for i = 1:size(myROIs,2)
            drawcircle('Position',myROIs{i},'Radius',radius);    
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

save('AllMyROIs','myROIs');
