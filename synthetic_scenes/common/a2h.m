% X = a2h(x [,w]) - affine to homogeneous coordinates
% 
% X = [x; 1 1 1 ...] 

% (c) T.Pajdla, www.neovision.cz, Sep 10, 2004
function X = a2h(x,w)

if nargin<2
    w = 1;
end
X = [x;w.*ones(1,size(x,2))];