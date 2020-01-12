%find Cdc42T peak and shift all profiles, so the peak is
%centered at (xnew,ynew) obtained as function parameters
% ee=1:1:N; 
% [X, Y]=meshgrid(ee);
% distm=sqrt((N/2 - X).^2 + (N/2 - Y).^2);
% indd=find((distm > dist0-1) & (distm < dist0+1));
% newi=randi([1, length(indd)]);
% [xnew0,ynew0] = ind2sub([N N],indd(newi));
x=importdata('input.dat'); 
xnew0=x(1);
ynew0=x(2);
% xnew0=floor(N*rand()); %chosen randomly
% ynew0=floor(N*rand()); %chosen randomly
%display(['Randomizing the initial peak position: (x0,y0) = (' num2str(xnew0) ',' num2str(ynew0) ')']);
display(['   Shifting the initial peak position: (x0,y0) = (' num2str(xnew0) ',' num2str(ynew0) ')']);
[indx,indy]=find(Cdc42T==max(max(Cdc42T)));
if (length(indx)>1) indx=indx(1); end
if (length(indy)>1) indy=indy(1); end
delx=(xnew0-indx);
dely=(ynew0-indy);
cdct=Cdc42T;cdcd=Cdc42D; gdicdc=GDI42;
bem=BemGEF; bem42=BemGEF42;vsn=vSNARE;reca=RecA; rec=Rec;newindx=0;newindy=0;
for i=1:N
    for j=1:N
        if (i+delx<0)
            newindx=N+i+delx;
        elseif (i+delx==0 | i+delx==N)
            newindx=N;
        else
            newindx=mod(i+delx,N);
        end
        if (j+dely<0)
            newindy=N+j+dely;
        elseif (j+dely==0 | j+dely==N)
            newindy=N;
        else
            newindy=mod(j+dely,N);
        end
        Cdc42T(  newindx,newindy)=cdct(  i,j);
        Cdc42D(  newindx,newindy)=cdcd(  i,j);
        GDI42(   newindx,newindy)=gdicdc(i,j);
        BemGEF(  newindx,newindy)=bem(   i,j);
        BemGEF42(newindx,newindy)=bem42( i,j);
        vSNARE(  newindx,newindy)=vsn(   i,j);
        RecA(    newindx,newindy)=reca(  i,j);
        Rec(     newindx,newindy)=rec(   i,j);
    end
end
display('   Shifted and PBC wrapped all concentration profiles...');
