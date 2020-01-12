function endo_cell = moving_patches_endo_radial(r_endo,X,Y,endo_cell,C)

% r_endo      = radious of exocytic patch
% [X,Y]       = centre of removal
% endo_cell   = endocytic patch data
% C           = circumference of cell

%% centre vesicle event
Xdx = (C/2) - X;
Ydy = (C/2) - Y;

Xc = C/2;
Yc = C/2;

x0 = mod(endo_cell.posx + Xdx,C);
y0 = mod(endo_cell.posy + Ydy,C);
ind = find(x0 == 0);
x0(ind) = C;
ind = find(y0 == 0);
y0(ind) = C;

% figure;hold;title('endo before centre')
% plot(X,Y,'go')
% plot(endo_cell.posx,endo_cell.posy,'o')
% figure;hold;title('endo after centre')
% plot(Xc,Yc,'go')
% plot(x0,y0,'ro')

x1 = NaN(1,endo_cell.num);
y1 = NaN(1,endo_cell.num);

% distance of endo patch from (Xc,Yc)
d0 = ( (Xc-x0).^2 + (Yc-y0).^2 ).^(0.5);
% distance of endo patch from (Xc,Yc) after movement
d1 = (d0.^2 - r_endo^2).^0.5;
d = abs(d1-d0);
d(find(d0>=C/2)) = 0;
d(find(d0<=r_endo)) = 0;

ind = find(Xc~=x0 & Yc ~= y0);

    y1p(ind) = y0(ind) + d(ind).*((1+(((Xc-x0(ind)).^2)./((Yc-y0(ind)).^2))).^(-0.5));
    y1m(ind) = y0(ind) - d(ind).*((1+(((Xc-x0(ind)).^2)./((Yc-y0(ind)).^2))).^(-0.5));

    x1p(ind) = x0(ind) - (y0(ind)-y1p(ind)).*(Xc-x0(ind))./(Yc-y0(ind));
    x1m(ind) = x0(ind) - (y0(ind)-y1m(ind)).*(Xc-x0(ind))./(Yc-y0(ind));
    
    for i = ind
        if (Xc-x1p(i))^2 + (Yc-y1p(i))^2 < (Xc-x1m(i))^2 + (Yc-y1m(i))^2
            x1(i) = x1p(i);
            y1(i) = y1p(i);
        else
            x1(i) = x1m(i);
            y1(i) = y1m(i);
        end
    end

 ind = setdiff([1:endo_cell.num],ind);

for i = ind
    
    if Xc == x0(i) & Yc == y0(i)
    
        y1(i) = y0(i);
        x1(i) = x0(i);
            
    elseif Xc == x0(i)
    
        x1(i) = x0(i);  

        y1p(i) = y0(i) + d(i);
        y1m(i) = y0(i) - d(i);

        if (Xc-x1(i))^2 + (Yc-y1p(i))^2 < (Xc-x1(i))^2 + (Yc-y1m(i))^2
            y1(i) = y1p(i);
        else
            y1(i) = y1m(i);
        end

    
    else  % Yc == y0
    
        y1(i) = y0(i);  

        x1p(i) = x0(i) + d(i);
        x1m(i) = x0(i) - d(i);

        if (Xc-x1p(i))^2 + (Yc-y1(i))^2 < (Xc-x1m(i))^2 + (Yc-y1(i))^2
            x1(i) = x1p(i);
        else
            x1(i) = x1m(i);
        end
    end
end

% figure;hold;title('endo')
% plot(Xc,Yc,'go')
% plot(x0,y0,'o')
% plot(x1,y1,'ro')
% for i = 1:endo_cell.num
%     plot([x0(i) x1(i)],[y0(i) y1(i)])
% end


% transform the moved patches back
endo_cell.posx = mod(x1 - Xdx,C);
endo_cell.posy = mod(y1 - Ydy,C);

end
