function [maxm,tmax] = get_degree(eqs);

for kk = 1:(length(eqs));
    mm = monomials(eqs(kk));
    maxm{kk} = max(mm')';
    if size(mm,1)>1
        tmax(kk) = max(sum(mm));
    else
        tmax(kk) = max(mm);
    end
end