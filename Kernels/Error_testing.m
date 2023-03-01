%Integral to check order of magnitude of the erros of the IB method
sig = tan(pi/6);
X_s = 1.0;
Y_s = sig*X_s;
l = 2;
h = 1/256;

f = @(x,y) 1/(4*h) .* (1 + cos(pi.*(x-X_s)/(2*h))).* 1/(4*h) .* (1 + cos(pi.*(y-Y_s)/(2*h))).*(y - sig*x)./(sqrt(sig*sig + 1));

ymin = @(x) sig*x;
xmax = X_s + l*h;
xmin = X_s - l*h;
ymax = sig*X_s + l*h;

integral2(f,xmin,xmax,ymin,ymax)
