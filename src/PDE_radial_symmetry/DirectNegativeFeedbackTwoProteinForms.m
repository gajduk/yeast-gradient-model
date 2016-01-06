classdef DirectNegativeFeedbackTwoProteinForms < ModelCore
    %DirectNegativeFeedbackTwoProteinForms
    %See https://www.overleaf.com/2061351vhqrxw#/5221863/ - label phosphatase negative feedback
    properties (Constant)
        R = 10^-5;%sphere radius
        D = 10*10^-13;%diffusion constant of the phosphorylated protein
        p0 = 300;%[nM]initial concentration of protein (all unphosphorylated)
        kp = 1;%s^-1 rate of de-phosphorylation
        ka = 1;%s^-1 rate of activation (phosphorylation)

        %phosphatase governed by Michaelis Menten kinetics
        k_cat = .05;%s^-1 catalytic activity of the phosphatases
        E0 = 75;%[nM] is the initial concentration of the phosphotases
        K_M = 100;%[nM] is the  Michaelis constant of the phosphatases
        A = 10;%strength of the feedback (a free parameter).
    end
    
    methods
       function [c,f,s] = pde_fun(self,x,t,u,DuDx)
            p = u(1);%unphosphorylated form of the protein
            pa = u(2);%active - phosphorylated form of the protein
            ea = self.ka*x^10/self.R^10;%effective rate of activation
            pr = self.phosphatase_rate(pa,self.k_cat,self.E0,self.K_M,self.A);%the activity of the phosphotases is proportionally to pa
            
            c = [1; 1];
            f = [self.D; self.D].*DuDx;
            s = [2*self.D/x; 2*self.D/x].*DuDx+[pr-ea*p;-pr+ea*p];
       end
       
       function res = phosphatase_rate(self,pa,k_cat,E0,K_M,A)
            %PHOSPHATASE_RATE - the activity of the phosphotases increases
            %proportionally with $u^p(r,t)$ where
            %k_cat = 1 [s^-1]  is the catalytic activity of the phosphatases, 
            %E_0 = 75 [nM] is the initial concentration of the phosphotases (note that the concentration doesn't change in the model), 
            %K_M = 100 [nM] is the  Michaelis constant and
            %A is the strength of the feedback (a free parameter).
            res = k_cat*((E0*pa)/(K_M+pa))*((1+A*pa*K_M)/(1+pa*K_M));
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

