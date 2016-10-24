clc;
clear;

ICPC = 0;
%   The Levenberg-Marquardt algorithm does not handle bound constraints
%   The trust-region-reflective algorithm requires at least as many equations as variables; aborting.
options = optimset('Display','iter','Algorithm','Levenberg-Marquardt');
x0 = [rand/4,rand/4,rand/4];
[x,fval,exitflag,output] = fsolve(@InformationContent,x0,options,ICPC);

x
[abs(x),1-sum(abs(x))]
fval
exitflag