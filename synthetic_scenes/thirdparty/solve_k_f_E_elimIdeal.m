% test 7pt k+f+E


function  [F,f,k] = solve_k_f_E_elimIdeal(x,xp)

% Viktor: Rescaling added
f0 = mean(sqrt(sum(x(:,1:2).^2,2)));
x(:,1:2) = x(:,1:2)/f0;

%fill the 7x12 matrix
% cplumns correspond to monomials
% f11, f12, f13, f21, f22, f23, f31, f32,  f33, f31*k, f32*k,f33*k,
M = [ x(:,1).*xp(:,1), x(:,1).*xp(:,2), x(:,1), x(:,2).*xp(:,1), x(:,2).*xp(:,2), x(:,2), xp(:,1), xp(:,2),  ones(size(x,1),1), xp(:,1).*(x(:,1).^2 + x(:,2).^2),xp(:,2).*(x(:,1).^2 + x(:,2).^2),x(:,1).^2 + x(:,2).^2];

%null space
n = null(M);

%GB solver
[a, b, c, d] = solver_E_f_l_sturmfels(n(:,1), n(:,2), n(:,3), n(:,4), n(:,5));

% reconstruct solutions
for i = 1:length(a)
    % 12dim vector
    FF = a(i)*n(:,1)+b(i)*n(:,2)+c(i)*n(:,3)+d(i)*n(:,4)+n(:,5);
    
    %fundamental matrix
    F{i} = reshape(FF(1:9)', 3,3) ;
    %radial distortion parameter
    k(i) = FF(10)/FF(7);
    
%     x11 = F{i}(1,1);
%     x12 = F{i}(1,2);
%     x13 = F{i}(1,3);
%     x21 = F{i}(2,1);
%     x22 = F{i}(2,2);
%     x23 = F{i}(2,3);
%     x31 = F{i}(3,1);
%     x32 = F{i}(3,2);
%     x33 = F{i}(3,3);
%     
    
    x11 = FF(1);
    x12 = FF(2);
    x13 = FF(3);
    x21 = FF(4);
    x22 = FF(5);
    x23 = FF(6);
    x31 = FF(7);
    x32 = FF(8);
    x33 = FF(9);
    
    f(i) = sqrt((x23*x31^2+x23*x32^2-2*x21*x31*x33-2*x22*x32*x33-x23*x33^2)/(2*x11*x13*x21+2*x12*x13*x22-x11^2*x23-x12^2*x23+x13^2*x23+x21^2*x23+x22^2*x23+x23^3));

    f(i) = f(i) * f0;
    k(i) = k(i) / f0^2;
    F{i} = F{i} * diag([f0 f0 1]);
end
