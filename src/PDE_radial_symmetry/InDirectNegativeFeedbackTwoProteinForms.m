classdef InDirectNegativeFeedbackTwoProteinForms < ModelCore
    %DirectNegativeFeedbackTwoProteinForms
    %See https://www.overleaf.com/2061351vhqrxw#/5221863/ - label phosphatase negative feedback
    properties (Constant)
        R = 10^-5;%sphere radius
        D = 5*10^-14;%diffusion constant of the phosphorylated protein
        p0 = 300;%[nM]initial concentration of protein (all unphosphorylated)
        phos0 = 50;%[nM] initial concentration of phospatases
        phos_max = InDirectNegativeFeedbackTwoProteinForms.phos0*4;%[nM] maximal concentration of phosphatases
        phos_regulation_rate = .002;% - the speed of phosphatase dynamics
        kp = .003;%s^-1 rate of de-phosphorylation
        ka = .008;%s^-1 rate of activation (phosphorylation)
        A = 1;% strength of the feedback (a free parameter).
    end
    
    methods
       function [c,f,s] = pde_fun(self,x,t,u,DuDx)

            p = u(1);%unphosphorylated form of the protein
            pa = u(2);%active - phosphorylated form of the protein
            phos = u(3);%phophatases
            ea = self.ka*x^10/self.R^10;%effective rate of activation
            pr = pa*self.kp*phos/self.phos0;

            c = [1; 1; 1];
            f = [self.D; self.D; self.D].*DuDx;
            s = [2*self.D/x; 2*self.D/x; 2*self.D/x].*DuDx;
            s = s+[ pr-ea*p;%non-active protein
                    -pr+ea*p;%active protein
                    self.phos_regulation_rate*(self.A*pa*(self.phos_max-phos)/(self.phos_max-self.phos0))];%phosphotases
       end
       
       function u0 = ic_fun(self,x)
            u0 = [self.p0; 0; self.phos0];
       end
       
       function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
            % pR is a constant and the initial concentration of phosphorylated protein at the membrane
            % D - diffusion constant
            % pl and ql correspond to the left boundary conditions (x = 0), and 
            % pr and qr correspond to the right boundary condition (x = R).
            pl = [0; 0; 0];
            ql = [1; 1; 1];
            pr = [0; 0; 0];
            qr = [1; 1; 1];
        end

    end
    
end

