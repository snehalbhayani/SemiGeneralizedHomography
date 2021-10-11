function sols = solver_radial_oneside_E_lambda_w(u1d,u2d);


r1d2 =  u1d(1,:).^2 + u1d(2,:).^2 ;

x = create_vars(10);
LAMBDA = x(9);
OMEGA = x(10);
F = reshape([x(1:8); 1],3,3);
U1 = [u1d(1:2,:);u1d(3,:)+LAMBDA*r1d2];
U2 = u2d;

for k = 1:7;
    eqs(k,1)=U1(:,k)'*F*U2(:,k);
end

%%
[C,monv] = polynomials2matrix(eqs);

sel1 = [3 4 6 7 9 10 2];
sel2 = [1 5 8 11 12];

BB = -C(:,sel1)\C(:,sel2);

y = create_vars(4);
% y(1) = x3 = f3
% y(2) = x6 = f6
% y(3) = x9 = lambda;
% y(4) = x10 = w;
monvy = [y(1)*y(3);y(1);y(2);y(3);1];
ff = BB*monvy;

F = [ff(1) ff(3) ff(5); ff(2) ff(4) ff(6);y(1) y(2) 1];
Kw = diag([1 1 y(4)]);
E = Kw*F;
tmp = E*E'; 
tmp2 = tmp(1,1)+tmp(2,2)+tmp(3,3);
eqs2 = 2*E*E'*E - tmp2*E;
eqs2 = [eqs2(:);det(E)];
eqs2(11) = y(2)*y(3)-ff(7);

% Divide equations 3 6 9 10 with w

for k = [3 6 9 10];
    oneeqs = eqs2(k);
    cc = coeffs(oneeqs);
    mm = monomials(oneeqs);
    mm(4,:) = mm(4,:)-1;
    eqs2(k) = multipol(cc,mm);
end

% substitute w^2 with w
for k = 1:11;
    oneeqs = eqs2(k);
    cc = coeffs(oneeqs);
    mm = monomials(oneeqs);
    mm(4,:) = mm(4,:)/2;
    eqs2(k) = multipol(cc,mm);
end

[cc,mm]=polynomials2matrix(eqs2);

%%

eqsn = eqs2;

[mm,nn] = get_degree(eqsn);
eqs3 = [];

N = 4;
M = 8;
NN = N*ones(4,1);
NN(1) = 3;
NN(2) = 4;
%NN(3) = 2;
%NN(4) = 4;
NN(3) = 4;
NN(4) = 2;
for i = 1:length(eqsn);
    
    
    
    multmons{i} = monvec(create_upto_degree(...
        max(NN-mm{i},0),...
        M-nn(i)));
    
    if 0
        if length(multmons{i}) == 42
            xx = [1:11,15,35:42];
        else
            xx = [1:31];
        end
        xx = [1:10,15:20,28:31];
        multmons{i} = multmons{i}(xx);
    end
    eqs3 = [eqs3;eqsn(i)*multmons{i}(:)];
    
end



settings.dim = 19;
settings.basis_size = 40;
settings.action_variable = 1;
[sols,stats] = polysolve(eqs3,settings);
