
function [ind,Xi,Yi,interpD] = sink_ind_02g (posx, posy, dx, N)

% number of points of endocytosis
nCell = length(posx);

ind     = NaN(nCell,4);
Xi      = NaN(2*nCell,2);
Yi      = NaN(2*nCell,2);
interpD = NaN(nCell,1);

for i=1:nCell       % for all the points of endocytosis    

    % it's possible that if the membrane gets small enough
    % that a vesicle position is outside N*dx
    % this is because we do not transform endo_cell_posx,posy
    % with the rest of the membrane
    % for now - when this happens wrap around
    if posx(i) > N*dx
        curr_t
        posx
        N*dx
        posx(i) = posx(i)-(N*dx);
    end
    if posy(i) > N*dx
        curr_t
        posy
        N*dx
        posy(i) = posy(i)-(N*dx);
    end    
        
    if (posx(i)/dx) - floor(posx(i)/dx) ~= 0        % not on an x grid line
        % finding vesicle's neighbouring x grid lines
        xm = floor(posx(i)/dx);             % xm = [0 --> N-1]
        xp = ceil(posx(i)/dx);              % xp = [1 --> N  ]  
        
        if (posy(i)/dx) - floor(posy(i)/dx) ~= 0    % not on a y grid line
            % finding vesicle's neighbouring y grid lines
            ym = floor(posy(i)/dx);         % ym = [0 --> N-1]
            yp = ceil(posy(i)/dx);          % yp = [1 --> N  ]  
            
            xVec = [xm,xp,xm,xp];
            yVec = [ym,ym,yp,yp];
            
            % positions of vesicle's neighbours, microns 
            % (don't wrap for interpolation)
            Xi((2*i)-1:2*i,:)  = dx*reshape(xVec,2,2)';
            Yi((2*i)-1:2*i,:)  = dx*reshape(yVec,2,2)';  
            interpD(i) = 2;
            
            % to define ind wrap around so all points in [1 --> N]            
            xVec(xVec==0) = N;
            yVec(yVec==0) = N;                        
            % index of vesicle's neighbours within distance dx
            ind(i,:) = sub2ind([N,N],xVec,yVec);    
            
        else                                        % on a y grid line
            y = posy(i)/dx;                 % y  = [1 --> N  ]
            
            xVec = [xm,xp];
            yVec = [y, y ];
            
            % positions of vesicle's neighbours, microns 
            % (don't wrap for interpolation) 
            Xi((2*i)-1,:) = dx*xVec;
            Yi((2*i)-1,:) = dx*yVec;    
            interpD(i)    = 1;            
            
            % to define ind wrap around so all points in [1 --> N]  
            xVec(xVec==0) = N;
            yVec(yVec==0) = N;       
            % index of vesicle's neighbours within distance dx            
            ind(i,1:2)    = sub2ind([N,N],xVec,yVec);
          
            
        end
    else                                            % on an x grid line
        if (posy(i)/dx) - floor(posy(i)/dx) ~= 0    % not on a y grid line
            ym = floor(posy(i)/dx);         % ym = [0 --> N-1]
            yp = ceil(posy(i)/dx);          % yp = [1 --> N  ] 
            x  = posx(i)/dx;                % x  = [1 --> N  ]
           
            xVec = [x ,x ];
            yVec = [ym,yp];

            % positions of vesicle's neighbours, microns 
            % (don't wrap for interpolation)             
            Xi((2*i)-1,:) = dx*xVec;
            Yi((2*i)-1,:) = dx*yVec;        
            interpD(i) = 1;
            
            % to define ind wrap around so all points in [1 --> N] 
            xVec(xVec==0) = N;
            yVec(yVec==0) = N;               
            ind(i,1:2)    = sub2ind([N,N],xVec,yVec);
      
        else                                        % on a y grid line
            ind(i,1) = sub2ind([N,N],posx(i)/dx,posy(i)/dx);
%             Xi  = NaN;
%             Yi  = NaN;               
            interpD(i) = 0;
        end
    end

    

end

return

