
%radius of vesicle
rv =r_endo;
% center of vesicle
Xc = C/2;
Yc = C/2;

% create regular grid we will finish with
x0 = [C/N:C/N:C];
y0 = [C/N:C/N:C];
x0 = repmat(x0,N,1);
y0 = repmat(y0',1,N);

x1 = NaN(N);
y1 = NaN(N);
x1p = NaN(N);
y1p = NaN(N);
x1m = NaN(N);
y1m = NaN(N);

% distance of grid point from (Xc,Yc)
d0 = ( (Xc-x0).^2 + (Yc-y0).^2 ).^(0.5);

% distance of irregular grid point from (Xc, Yc)
d1 = real(( (d0.^2) + (rv^2) ).^(0.5));
d = abs(d0-d1);
d(find(d0>=C/2)) = 0;
% figure;surf(d)

ind = find(Xc~=x0 & Yc ~= y0);

    y1p(ind) = y0(ind) + d(ind).*((1+(((Xc-x0(ind)).^2)./((Yc-y0(ind)).^2))).^(-0.5));
    y1m(ind) = y0(ind) - d(ind).*((1+(((Xc-x0(ind)).^2)./((Yc-y0(ind)).^2))).^(-0.5));

    x1p(ind) = x0(ind) - (y0(ind)-y1p(ind)).*(Xc-x0(ind))./(Yc-y0(ind));
    x1m(ind) = x0(ind) - (y0(ind)-y1m(ind)).*(Xc-x0(ind))./(Yc-y0(ind));
    
    for j = 1:length(ind)
        i = ind(j);
        if (Xc-x1p(i))^2 + (Yc-y1p(i))^2 > (Xc-x1m(i))^2 + (Yc-y1m(i))^2
            x1(i) = x1p(i);
            y1(i) = y1p(i);
        else
            x1(i) = x1m(i);
            y1(i) = y1m(i);
        end
    end

ind = setdiff([1:N*N],ind);

for j = 1:length(ind)
    i = ind(j);    
    if Xc == x0(i) & Yc == y0(i)
    
        y1(i) = y0(i);
        x1(i) = x0(i);
            
    elseif Xc == x0(i)
    
        x1(i) = x0(i);  

        y1p(i) = y0(i) + d(i);
        y1m(i) = y0(i) - d(i);

        if (Xc-x1(i))^2 + (Yc-y1p(i))^2 > (Xc-x1(i))^2 + (Yc-y1m(i))^2
            y1(i) = y1p(i);
        else
            y1(i) = y1m(i);
        end

    
    else  % Yc == y0
    
        y1(i) = y0(i);  

        x1p(i) = x0(i) + d(i);
        x1m(i) = x0(i) - d(i);

        if (Xc-x1p(i))^2 + (Yc-y1(i))^2 > (Xc-x1m(i))^2 + (Yc-y1(i))^2
            x1(i) = x1p(i);
        else
            x1(i) = x1m(i);
        end
    end
end

% figure;hold
% for i = 1:N
%     for j = 1:N        
%         plot([x0(i,j) x1(i,j)],[y0(i,j) y1(i,j)])
%     end
% end
% plot(x0(1:N,1:N),y0(1:N,1:N),'r.')
% axis([45 55 45 55])


x0 = [0:C/N:C];
y0 = [0:C/N:C];
% create regular grid
x0 = repmat(x0,N+1,1);
y0 = repmat(y0',1,N+1);





