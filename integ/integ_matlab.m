clear, close all; clc;

syms x t
b = 0; a = -inf;
f(x) = exp(-x^2);

res1 = int(f,a,b);

tb = 0; ta = -1;
ft(t) = f(b-t/(1+t))*1/(1+t)^2;

res2 = int(ft,ta,tb);

disp('-----------------------------------');
disp('integral  of'); f
disp([' from ',num2str(a),' to ',num2str(b)]); 
disp(res1);
disp('-----------------------------------');
disp('integral  of'); ft
disp([' from ',num2str(ta),' to ',num2str(tb)]); 
disp(res2);