%Solve_2f computes focal lengths from a fundamental matrix.
%
%[f1,f2,C,omegas] = Solve_2f(F)
%
%   Given the fundamental matrix, F, that captures the epipolar geometry of
%   two cameras, this function returns the two focal length of the two
%   cameras.
%
%   Although the fundamental matrix embodies both the intrinsic parameters
%   (cameras' focal lengths and principal points) and extrinsic parameters
%   (relative orientation and translation between the two cameras), it only
%   has 7 degrees of freedom.  This limits the number of recoverable
%   parameters to 7.  There are 5 extrinsic parameters (3 rotation angles
%   and 2 components of the translation vector - note that the translation
%   vector is recoverable only up to an unknown scale).  So only 2
%   intrinsic parameters can be estimated.
%
%   The code implemented here is based on the algorithm described in the
%   following paper:
%   G. N. Newsam, D. Q. Huynh, M. J. Brooks, H.-P. Pan.
%   "Recovering Unknown Focal Lengths in Self-Calibration: An Essentially
%   Linear Algorithm and Degenerate Configurations."
%   International Archives of Photogrammetry and Remote Sensing, 
%   vol XXXI, part B3, commission III, pp. 575-580, Vienna, Austria,
%   9-19 July 1996.
%   
%   The algorithm assumes that the principal points of both cameras are
%   known and that the origin of the image coordinate system has been fixed
%   at the locations of the known principal points when the fundamental
%   matrix, F, was computed.
%
%   If the input F matrix was poorly estimated then it is possible that
%   the function failed to compute the focal length.  In this case, a
%   message would be displayed and the focal length would be set to 0.
%
%   Note that if the two cameras' optical axes and the baseline are
%   coplanar in space then it would not be possible to recover the two
%   unknown focal lengths.  However, if the two focal lengths are identical
%   (you would then use the function Solve_f instead) then it would be
%   possible to recover the single focal length shared by the two cameras.
%
%   Use Solve_f.m if the two cameras have the same focal length.
%
%   Input parameter:
%   F     the input fundamental matrix
%
%   Output parameters:
%   f1, f2:    the two estimated focal lengths.
%   C:         the matrix as described in Eqs(26-27).
%   sigmas:    an array containing the 3 singular values of F.
%
%SEE ALSO Solve_f
%
%Du Huynh, created as a Maple file on 31 March 1996.
%Du Huynh, converted to Matlab in August 1999.

function [f1,f2,C,sigmas] = Solve_2f(F)

% normalize the given fundamental matrix
Fhat = normF(F);

% compute the SVD of Fhat
[U,S,V] = svd(Fhat);
sigma1 = S(1,1);
sigma2 = S(2,2);
sigma3 = S(3,3);
   
% u1, u2, u3 are , respectively, the 1st, 2nd, 3rd columns of U.
% Similarly for V and F
u1 = U(:,1);   u2 = U(:,2);   u3 = U(:,3);
v1 = V(:,1);   v2 = V(:,2);   v3 = V(:,3);
f1 = Fhat(:,1);   f2 = Fhat(:,2);   f3 = Fhat(:,3);
i3 = [0;0;1];
% construct the C matrix as described in Eqs (26-27)
C0 = makeC(u1, u2, u3, f1, f2, f3, i3);

% construct vector s that contains the sigmas
s = [ sigma1^2; 0; sigma2^2 ];

% check det(C)
if cond(C0) < 1e10
	omegas = C0 \ s;

	mu = -omegas(1);
	lambda = omegas(3);
	nu = omegas(2) / lambda;

	f1 = sqrt( 1/(1 + mu) );
   f2 = sqrt(1 + nu);
   if abs(imag(f1)) > 1e-5
    %  fprintf('Solve_2f: estimated f1 has an imaginary component\n');
   end
   if abs(imag(f2)) > 1e-5
   %   fprintf('Solve_2f: estimated f2 has an imaginary component\n');
   end
else
   %fprintf('Solve_2f: matrix C is singular.  ');
   %fprintf('Stereo configuration is degenerate!\n');
	f1 = 0;
	f2 = 0;
end

if nargout > 3
   C = C0;
   sigmas = diag(S);
elseif nargout > 2
   C = C0;
end
return

% ------------------------------------------------------------------

% function C = makeC(u1,u2,u3,f1,f2,f3,i3)
%
% This function is called by Solve_2f.m only  and is not expected
% to be needed by any other Matlab routines.
% Arguments:
%  u1,u2,u3: the 1st, 2nd, 3rd columns of the U matrix where U is obtained
%            from the singular value decomposition of the fundamental matrix.
%  f1,f2,f3: the 1st, 2nd, 3rd columns of the fundamental matrix
%  i3:       this is simply [0,0,1]'
% Given the above arguments, the function returns the C matrix as described
% in Equations (26-27) of the Newsam et al's ISPRS 96 paper.
%
% Du Huynh, created on Sun, 30 March 1996 as a Maple program
% Du Huynh, converted in August 1999 to Matlab

function C = makeC(u1, u2, u3, f1, f2, f3, i3)

% Note that I have not attempted to optimise the code below.
C = zeros(3,3);
C(1,:) = [dot(u1,f3)^2 dot(u1,i3)^2+dot(u3,i3)^2 1];
C(2,:) = [dot(u2,f3)*dot(u1,f3) dot(u1,i3)*dot(u2,i3) 0];
C(3,:) = [dot(u2,f3)^2 dot(u2,i3)^2+dot(u3,i3)^2 1];
return
