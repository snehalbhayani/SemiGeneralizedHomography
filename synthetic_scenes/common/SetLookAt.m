% Matlab code by Zuzana Kukelova 2006-01-10
% Generate camera matrix
% From MartinB
% Input:    pos - Position of Camera
%           at -  Look at - point on which is camera looking at
%           up -  Up vector
% Output:   Camera = normalized camera matrix K^-1*P = [R | t]


function P = SetLookAt(pos, at, up)

  %direction vector - z
  n = at - pos;  
  n = n./(sqrt(sum(n.*n)));
  
  %right vector - x
  u = cross(n,up);	
  u = u./(sqrt(sum(u.*u)));
  
  % up vector - y
  v = cross(u,n);	
  v = -v./(sqrt(sum(v.*v)));
  
  t = [ -dot(u,pos),-dot(v,pos),-dot(n,pos)];

  P = [u;v;n;t];
  P = P';
