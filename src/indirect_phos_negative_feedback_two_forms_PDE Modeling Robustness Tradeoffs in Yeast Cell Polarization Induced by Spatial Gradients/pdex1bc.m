function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
% pR is a constant and the initial concentration of phosphorylated protein at the membrane

% pl and ql correspond to the left boundary conditions (x = 0), and 
% pr and qr correspond to the right boundary condition (x = R).

pl = [0; 0; 0];
ql = [1; 1; 1];

pr = [0; 0; 0];
qr = [1; 1; 1];


