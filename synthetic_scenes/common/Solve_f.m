%Solve_f computes focal length from a fundamental matrix.
%
%f = Solve_f(F)
%   Given the fundamental matrix, F, that captures the epipolar geometry of
%   two cameras, both of which have the same but unknown focal length, this
%   function returns the focal length of the two cameras.
%
%   Although the fundamental matrix embodies both the intrinsic parameters
%   (cameras' focal lengths and principal points) and extrinsic parameters
%   (relative orientation and translation between the two cameras), it only
%   has 7 degrees of freedom.  This limits the number of recoverable
%   parameters to 7.  There are 5 extrinsic parameters (3 rotation angles
%   and 2 components of the translation vector - note that the translation
%   vector is recoverable only up to an unknown scale).  So only 2
%   intrinsic parameters can be estimated, namely the cameras' focal
%   lengths.
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
%   Use Solve_2f.m if the two cameras are not identical.
%
%SEE ALSO Solve_2f
%
%Du Huynh, created as a Maple file on 30 March 1996.
%Du Huynh, converted to Matlab in August 1999.

function f = Solve_f(F)

% normalize the given fundamental matrix
Fhat = normF(F);

% compute the SVD of Fhat
[U,S,V] = svd(Fhat);
sigma1 = S(1,1);
sigma2 = S(2,2);
sigma3 = S(3,3);
   
% u1, u2, u3 are, respectively, the 1st, 2nd, 3rd columns of U.
% Similarly for V and F
u1 = U(:,1);   u2 = U(:,2);   u3 = U(:,3);
v1 = V(:,1);   v2 = V(:,2);   v3 = V(:,3);
f1 = Fhat(:,1);   f2 = Fhat(:,2);   f3 = Fhat(:,3);
i3 = [0;0;1];

% construct the quadratic in (39):
% coeff0 + coeff1*mu + coeff2*F33*mu^2 = 0
% we need to solve for mu
coeff0 = sigma1^2 - sigma2^2;
coeff1 = ( dot(u1,i3)^2 + dot(v1,i3)^2 ) * sigma1^2 - ...
   ( dot(u2,i3)^2 + dot(v2,i3)^2 ) * sigma2^2;
coeff2 = ( dot(u1,i3)*dot(v1,i3)*sigma1 - dot(u2,i3)*dot(v2,i3)*sigma2 ) * ...
   Fhat(3,3);

% solve the quadratic equation for mu.  The roots are either both real or
% both complex.
mu = roots([coeff2, coeff1, coeff0]);

% Ideally, we should have two distinct real roots for mu: mu1 and mu2 such
% that the inequality condition  mu1 < -1 < mu2  holds.  However, the
% input F matrix could be poorly estimated and the focal length cannot be
% estimated.


if abs(imag(mu(1))) < 1e-7
   % Both roots are real (imaginary component is negligible)
   if abs(real(mu(1)-mu(2))) < 1e-5
      % roots are identical
   %   fprintf('Solve_f: identical roots for Eq (39) (see paper)\n');
      f = 0;
      if abs(imag(f)) > 1e-5
         % computed focal length is complex because mu is small, negative
   %      fprintf('Solve_f: focal length was complex!\n');
         f = 0;
      end
   else
      % we have two distinct real roots
      est_fs = sqrt( 1./(mu + 1) );

      % one of these two computed focal lengths must be complex (pure
      % imaginary) while the other is real.  We take the real root as
      % the answer.
      if all(imag(est_fs) ~= 0)
         % both are real
    %     fprintf('Solve_f: found 2 distinct focal lengths!\n');
         f = 0;
      elseif all(imag(est_fs) == 0)
    %     fprintf('Solve_f: found 2 pure imaginary focal lengths!\n');
         f = 0;
      else 
         f = est_fs(find(imag(est_fs) == 0));
      end
   end
else
 %  fprintf('Solve_f: mu, the roots to Eq(39) are complex (conjugate pair)!\n');
   f = 0;
end

return
