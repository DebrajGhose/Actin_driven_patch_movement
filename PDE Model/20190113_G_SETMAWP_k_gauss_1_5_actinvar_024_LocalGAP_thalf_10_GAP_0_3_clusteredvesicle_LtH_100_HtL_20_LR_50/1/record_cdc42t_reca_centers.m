reci=Cdc42T+BemGEF42;
[xr,yr]=find(reci>0.9*max(max(reci)));
indr=find(reci>0.9*max(max(reci)));
%[xr,yr]=find(reci==max(max(reci)));
if((max(xr)- min(xr))>N/2) %space between found patches is larger than half of the box->the patch is going over the boundary
    is1=find(xr<N/2); % bins found on the left side
    is2=find(xr>=N/2); % bins found on the right side
    if (length(is2)>length(is1)) %majority of bins are on the right, easier to shift to the right
        xmean=(sum(reci(sub2ind(size(reci),xr(is1),yr(is1))).*(xr(is1)+N))+sum(reci(sub2ind(size(reci),xr(is2),yr(is2))).*xr(is2)))/sum(reci(indr));
        if (xmean>N) xmean=xmean-N; end %the larger part of the patch is on the left.
    else
        xmean=(sum(reci(sub2ind(size(reci),xr(is1),yr(is1))).*xr(is1))+sum(reci(sub2ind(size(reci),xr(is2),yr(is2))).*(xr(is2)-N)))/sum(reci(indr));
        if (xmean<0) xmean=xmean+N; end %on the right. Note: we have to pick one side anyway, including the case s1=s2.
    end
else xmean=sum(reci(indr).*xr)/sum(reci(indr)); %the patch is grouped in one part of the box, use the geometric average.
end
if((max(yr)- min(yr))>N/2) %space between found patches is larger than half of the box->the patch is going over the boundary
    is1=find(yr<N/2); % bins found on the top
    is2=find(yr>=N/2); % bins found on the bottom
    if (length(is2)>length(is1))
        ymean=(sum(reci(sub2ind(size(reci),xr(is1),yr(is1))).*(yr(is1)+N))+sum(reci(sub2ind(size(reci),xr(is2),yr(is2))).*yr(is2)))/sum(reci(indr));
        if (ymean>N) ymean=ymean-N; end%the larger part of the patch is on the left.
    else
        ymean=(sum(reci(sub2ind(size(reci),xr(is1),yr(is1))).*yr(is1))+sum(reci(sub2ind(size(reci),xr(is2),yr(is2))).*(yr(is2)-N)))/sum(reci(indr));
        if (ymean<0) ymean=ymean+N; end%on the right. Note: we have to pick one side anyway, including the case s1=s2.
    end
else ymean=sum(reci(indr).*yr)/sum(reci(indr)); end
c42cen=[c42cen; [xmean,ymean]];


%% collect positions of actin cables

if ~exist('storecables','var')
    storecables = [];
end

if numel(cable)>0
    storecables =  [storecables ;    padarray(cable,[0,20-numel(cable)],-1,'post') ]; % when you don't have all ten cables, just pad with -1
else
    storecables = -1*ones(1,20); %when no cables are present, just assign -1
end

%% collect timepoints

if ~exist('timepoints','var')
    timepoints = [];
end

timepoints = [timepoints,curr_t]; % collect timepoints