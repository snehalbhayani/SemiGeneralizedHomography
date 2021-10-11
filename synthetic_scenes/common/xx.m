function s = xx(x,s)
% 
% x ...  3 x n matrix [a b c; ...]
% s ... {[0 -a b; a 0 -c; b c 0] ... }
if nargin>1
    s = [    0    -x(3)  x(2)
            x(3)    0    -x(1)
           -x(2)  x(1)    0 ];   
elseif size(x,2)==1
    s = [    0    -x(3)  x(2)
           x(3)    0    -x(1)
          -x(2)  x(1)    0 ];   
else
    for i=1:size(x,2)
        s{i} = [    0    -x(3,i)  x(2,i)
                  x(3,i)    0    -x(1,i)
                 -x(2,i)  x(1,i)    0 ];   
    end
end