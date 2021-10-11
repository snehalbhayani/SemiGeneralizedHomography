function [r q] = rq(A)

[q r] = qr(A(end:-1:1,end:-1:1)');
q = q(end:-1:1,end:-1:1)';
r = r(end:-1:1,end:-1:1)';

s = diag(sign(diag(r)));
r = r*s;
q = s*q;

if det(q)<0
	r = -r;
	q = -q;
end