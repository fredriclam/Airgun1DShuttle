function [Y] = LeightonsRK4(F,T,k,y0)

N = ceil(T/k);
t = 0;
y = y0;
Y = zeros(length(y),N);
for i = 1:N
	Y(:,i) = y;
	k1 = F(t,y);
	k2 = F(t+k/2,y+ k*k1/2);
	k3 = F(t+k/2,y+ k*k2/2);
	k4 = F(t+k, y + k*k3);
	y = y + k/6*(k1 + 2*k2 + 2*k3 + k4);
	t = t+k;
end
Y(:,i) = y;
