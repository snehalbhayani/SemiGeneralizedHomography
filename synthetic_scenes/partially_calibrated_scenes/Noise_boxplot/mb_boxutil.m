function [outlier,loadj,upadj,yy,q1,q3,n2,med,n1] = mb_boxutil(x,notch,whis,whissw)
%BOXUTIL Produces a single box plot.

% define the median and the quantiles
pctiles = prctile(x,[25;50;75]);
q1 = pctiles(1,:);
med = pctiles(2,:);
q3 = pctiles(3,:);

% find the extreme values (to determine where whiskers appear)
vhi = q3+whis*(q3-q1);
upadj = max(x(x<=vhi));
if (isempty(upadj)), upadj = q3; end

vlo = q1-whis*(q3-q1);
loadj = min(x(x>=vlo));
if (isempty(loadj)), loadj = q1; end

outlier = x<loadj | x > upadj;
yy = x(outlier);

if whissw == 0
   upadj = max(upadj,q3);
   loadj = min(loadj,q1);
end

if notch
    n1 = med + 1.57*(q3-q1)/sqrt(length(x));
    n2 = med - 1.57*(q3-q1)/sqrt(length(x));
    %prevent notches from extending past edge of box
    if n1>q3, n1 = q3; end
    if n2<q1, n2 = q1; end
else
    n1=med;
    n2=med;
end

end