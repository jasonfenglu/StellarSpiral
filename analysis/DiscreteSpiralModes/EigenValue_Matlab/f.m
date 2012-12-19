function f = f(x)
i = length(x);
f = zeros(size(x));

for j = 1 : i
	if x(j)<1.
		f(j) = 1.-(1.+4.*x(j)).*(1.-x(j))^4.;
	else
		f(j) = 1.;
end
end
