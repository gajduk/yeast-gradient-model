classdef BasicModelTwoProteinForms < ModelCore
    %BasicModelTwoProteinForms - See https://www.overleaf.com/2061351vhqrxw#/5221863/ - label initial
    
    properties (Constant)
        R = 10^-5;%sphere radius
        D = 5*10^-13;%diffusion constant of the phosphorylated protein
        p0 = 100;%[nM[ initial concentration of protein (all unphosphorylated)
        kp = 1;%[s^-1] rate of de-phosphorylation
        ka = 1;%[s^-1] rate of activation (phosphorylation)
    end
    
    methods
       function [c,f,s] = pde_fun(self,x,t,u,DuDx)
            p = u(1);%unphosphorylated form of the protein
            pa = u(2);%active - phosphorylated form of the protein
            ea = self.ka*self.get_kinease_activity_normalizer(x,self.R);%effective rate of activation
            c = [1; 1];
            f = [self.D; self.D].*DuDx;
            s = [2*self.D/x; 2*self.D/x].*DuDx+[self.kp*pa-ea*p;-self.kp*pa+ea*p];
       end
       
       function u0 = ic_fun(self,x)
            u0 = [self.p0; 0];
       end
       
       function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
            % pR is a constant and the initial concentration of phosphorylated protein at the membrane
            % D - diffusion constant
            % pl and ql correspond to the left boundary conditions (x = 0), and 
            % pr and qr correspond to the right boundary condition (x = R).
            pl = [0; 0];
            ql = [1; 1];
            pr = [0; 0];
            qr = [1; 1];
        end

    end
    
end

