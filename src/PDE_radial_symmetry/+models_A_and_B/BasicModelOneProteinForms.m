classdef BasicModelOneProteinForms < ModelCore
    %BasicModelOneProteinForms - See Modeling Robustness Tradeoffs in Yeast Cell Polarization Induced by Spatial Gradients
    %This is the model described in eq. 1
    
    properties (Constant)
        R = 10^-5;%sphere radius
        D = 1*10^-13;%diffusion constant of the phosphorylated protein
        pR = 1;%initial condition of the phosphorylated protein at the membrane p(x=R,t=0)
        kp = 0;%rate of de-phosphorylation
    end
    
    methods
       function [c,f,s] = pde_fun(self,x,t,u,DuDx)
          c = 1;
          f = self.D*DuDx;
          s = 2*self.D*DuDx/x-self.kp*x;
       end
       
       function u0 = ic_fun(self,x)
            u0 = x*.1;
       end
       
       function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
            % pl and ql correspond to the left boundary conditions (x = 0), and 
            % pr and qr correspond to the right boundary condition (x = R).
            pl = 0;
            ql = 1;
            pr = 0;
            qr = 1;
       end

    end
    
end

