function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);
