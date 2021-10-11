function F = normF(F)

norm_f3 = norm(F(:,3), 2);
   
% if the norm of f3 is too small for safe scaling, then do not scale
if norm_f3 > 10^(-5)
   scaledF = F / norm_f3;
else
   scaledF = F;
end

