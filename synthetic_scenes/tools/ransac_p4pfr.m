function [ f,R,t,k,best_inl ] = ransac_p4pfr( u, U, tol, iter )

f = [];
R = [];
t = [];
k = [];
best_inl = 0;

d = sum(u(1:2,:).^2);

for iter = 1:iter
    sample = randperm(size(U,2),4);
    
    warning('off','MATLAB:rankDeficientMatrix');
    [ff,Rt,kk] = solver_p4p_fr_null_40x50(u(:,sample),U(:,sample));
    warning('on','MATLAB:rankDeficientMatrix');
    
    for i = 1:length(ff)
        
        P = diag([ff(i) ff(i) 1])*Rt(:,:,i);
        
        % Compute projection
        proj = P*U;
        proj = proj(1:2,:) ./ proj([3;3],:);
        proj = proj .* (1+kk(i)*[d;d]);
        
        err = sqrt(sum((u-proj).^2));
        inl = nnz(err < tol);
    
        if best_inl < inl
            best_inl = inl;
            f = ff(i);
            R = Rt(:,1:3,i);
            t = Rt(:,4,i);
            k = kk(i);
        end
    end
end

end

