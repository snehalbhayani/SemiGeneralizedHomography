% create a set of k different random numbers from interval 1..n
% 

function [r] = randnk(n, k)

    if k > n
        
        k = n;
        r = randperm(n);
        return;
    end
    
    r = zeros(1, k);
    i = 1;
    while i <= k
        
        val = ceil(rand()*n);
        if sum(r(1:i)==val) == 0
            
            r(i) = val;
            i = i+1;
        end
    end
end