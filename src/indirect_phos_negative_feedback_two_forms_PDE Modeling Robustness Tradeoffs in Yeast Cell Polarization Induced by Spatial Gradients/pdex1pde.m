function [c,f,s] = pdex1pde(x,t,u,DuDx,D,kp,ka,R,A,phos_max,phos_regulation_rate,phos0)
%x - radius, spatial coordinate of derivation
%u - a column vector of the unphosphorylated and the phosphorylated form of
%a protein integral u dx = const. because the proeins can only pass from one form to another  
%D - diffusion constant, assumed same for both forms of the protein
%kp - basic rate of phosrylation in the cytosol
%ka - rate of activation which hapens close to the cytosol i.e. effective
%rate of activation = ka*(R-r)**2
%R total cell radius
%A - strength of the negative feedback
%phos_max - maximal concentration of phosphatases
%phos_regulation_rate - the speed of phosphatse dynamics
global pr_list
global t_list

p = u(1);%unphosphorylated form of the protein
pa = u(2);%active - phosphorylated form of the protein
phos = u(3);%phophatases
ea = ka*x^10/R^10;%effective rate of activation
pr = pa*kp*phos/phos0;


if x > 0.98*R
    pr_list = [pr_list,phos];
    t_list = [t_list,t];
end
    
c = [1; 1; 1];
f = [D; D; D].*DuDx;
s = [2*D/x; 2*D/x; 2*D/x].*DuDx+[pr-ea*p;-pr+ea*p;phos_regulation_rate*((1+A*pa)*(phos_max-phos)/(phos_max-phos0)-1)];
