% select 4 non-coplanar points among a set of all points
%
% by Martin Bujnak (c) 2010

function [x X e] = Select5NoCoplanar(m, M, eps)

    % select 
    watchDog = 0;
    cnt = size(m, 2);
    Xb = [];
    Xe = 0;
    while (watchDog < 1000)

        s = randnk(cnt, 5);
        
        x = m(:,s);
        X = M(:,s);

        v1 = X(:,3)-X(:,1);
        v2 = X(:,2)-X(:,1);
        
        d1 = norm(v1);
        d2 = norm(v2);
        
        n = cross(v1(1:3),v2(1:3));
        n = n./norm(n);
        d = -dot(n,X(1:3,1));
        dist = abs(dot(n,X(1:3,4)) + d);
        dist2 = abs(dot(n,X(1:3,5)) + d);

        dd = max(d1, d2);
        if (dist > dd*eps)  && (dist2 > dd*eps)

            X(4,:) = 1;
            x(3,:) = 1;
            e = 0;
            return;
        end
        
        if (dd > 0) && (dist/dd > Xe) && (dist2/dd > Xe)
            
            Xe = dd/dist;
            Xb = s;
        end
            
        watchDog = watchDog + 1;
    end

    X = M(:, Xb);
    x = m(:, Xb);
    X(4,:) = 1;
    x(3,:) = 1;

    e = 1;
end