function f_w1=fw01(sw,A,muw)
global swi sor muo
syms  swd sw1 
swd=(sw1-swi)/(1-swi-sor);
kro=1.*(1-swd).^2;krw=1.*(swd.^2);
f_w=1/(1+kro.*muw/(muo.*krw));
f_w2=diff(f_w,sw1)-f_w/(sw1-A); 
f_w1=subs(f_w2,sw);
