function [c,f,s] = pdex1pde(x,t,u,DuDx,D,kp,ka,R)
%x - radius, spatial coordinate of derivation
%u - a column vector of the unphosphorylated and the phosphorylated form of
%a protein integral u dx = const. because the proeins can only pass from one form to another  
%D - diffusion constant, assumed same for both forms of the protein
%kp - rate of phosrylation in the cytosol, the rate is constant because we
%assume the phosphatases are spread unfiormly in the cytosol
%ka - rate of activation which hapens close to the cytosol i.e. effective
%rate of activation = ka*(R-r)**2
%R total cell radius

p = u(1);%unphosphorylated form of the protein
pa = u(2);%active - phosphorylated form of the protein
ea = ka*x^10/R^10;%effective rate of activation
c = [1; 1];
f = [D; D].*DuDx;
s = [2*D/x; 2*D/x].*DuDx+[kp*pa-ea*p;-kp*pa+ea*p];
