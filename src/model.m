function [dy] = model(t,y,p,D,k,n)
%A simple n-compartment model
%The first n variables are concentrations of UNphophorylated protein `Fus3`
%The last n variables are concentrations of phophorylated protein `Fus3_p`
%The first compartment is the schmoo tip where phosphorylation happens with
%rate `k` and no DEphosphorylation
%In the rest of the compartments there is no phophorylation only
%DEphosphorylation with rate `p`
%Each compartment is subject to diffusion with neighbouring compartments

dy = zeros(2*n,1);%the change in protein concentration

%the intermediate compartments
for i=2:n-1
    dy(i) = D*(y(i-1)+y(i+1)-2*y(i))+p*y(i+n);%UNphosphorylated Fus3
    dy(n+i) = D*(y(n+i-1)+y(n+i+1)-2*y(n+i))-p*y(i+n);%phosphorylated Fus3_p
end
%The schmoo tip
dy(1) = D*(y(2)-y(1))-k*y(1);
dy(n+1) = D*(y(n+2)-y(n+1))+k*y(1);
%The last compartment
dy(n) = D*(y(n-1)-y(n))+p*y(n+n);
dy(n+n) = D*(y(n+n-1)-y(n+n))-p*y(n+n);
end

