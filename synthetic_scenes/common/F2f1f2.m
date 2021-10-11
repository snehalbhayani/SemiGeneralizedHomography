% f = F2f1f2(F,[p,mth]) - Focal lengths from Fundamental matrix
%
% F   = Fundamental matrix with diagonal K1 = diag([f1 f2 1]), K2 = diag([f2 f2 1])
% p   = {p1 p2} princial points ({[0;0] [0;0]} implicit)
% mth = method
%      'Bougnoux' - Bougnoux formula from HZ (implicit)
%
% T. Pajdla, pajdla@cvut.cz, 2016-08-28
function f = F2f1f2(F,p,mth)
if nargin>0
    if nargin<3
        mth = 'Bougnoux';
    end
    if nargin<2
        p = {[0;0] [0;0]};
    end
    switch mth
        case 'Bougnoux'
            % epipoles
            [e2,~,e1] = svd(F);
            e1 = e1(:,3);
            e2 = e2(:,3);
            II = diag([1 1 0]);
            p1 = a2h(p{1});
            p2 = a2h(p{2});
            f(1) = sqrt((-p2'*xx(e2)*II*F *(p1*p1')*F'*p2)/(p2'*xx(e2)*II*F*II*F'*p2));
            f(2) = sqrt((-p1'*xx(e1)*II*F'*(p2*p2')*F *p1)/(p1'*xx(e1)*II*F'*II*F *p1));
    end
else % unit tests
    p = {[500;600] [500;400]};
    fg = [1100 900];
    K1 = [fg(1)  0     p{1}(1)
          0      fg(1) p{1}(2)
          0      0     1];
    K2 = [fg(2) 0     p{2}(1)
          0     fg(2) p{2}(2)
          0     0     1];
    % Rotations
    R1 = eye(3);
    R2 = a2r(rand(3,1),pi/10);
    % Projection matrices
    P1 = K1*R1*[eye(3)       [-1000;0;0]];
    P2 = K2*R2*[eye(3)  1000*[rand(2,1);0]];
    F = PP2F(P1,P2);
    f = F2f1f2(F,p);
    f = all(abs(f-fg)<1e-8);
end