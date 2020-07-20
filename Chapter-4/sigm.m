function y = sigm(v)

v0=6;
e0=2.5;
r=0.56;

y = (2*e0)./(1+exp(r*(v0-v)));

end