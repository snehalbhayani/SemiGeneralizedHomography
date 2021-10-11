function xf = pflat(x)
xf = x ./ (ones(size(x,1),1)*x(end,:));