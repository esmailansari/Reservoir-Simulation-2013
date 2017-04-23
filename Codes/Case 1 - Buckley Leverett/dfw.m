function f_w=dfw(sw,muw)
global swi sor muo
syms  sw1 
swd=(sw1-swi)/(1-swi-sor);
kro=1*(1-swd).^2;krw=1*(swd.^2);
f_w0=1/(1+kro*muw/(muo*krw));
fwhelp=diff(f_w0,sw1);
f_w=subs(fwhelp,sw);