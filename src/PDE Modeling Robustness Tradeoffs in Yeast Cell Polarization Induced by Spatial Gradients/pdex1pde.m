function [c,f,s] = pdex1pde(x,t,u,DuDx,D,kp)
%D - diffusion constant
%kp rate of phosrylation in the cytosol, the rate is constant because we
%assume the phosphatases are spread unfiormly in the cytosol
c = 1;
f = D*DuDx;
s = 2*D*DuDx/x-kp*x;