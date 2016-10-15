classdef MMBasicModelTwoProteinForms < ModelCore
    %BasicModelTwoProteinForms
    %with michaelis Menten kinetics
    %See https://www.overleaf.com/2061351vhqrxw#/5221863/ - label michaelis
    %menten first TAC
    
    properties
        R = 10^-5;%sphere radius
        D = 5*10^-13;%diffusion constant of the phosphorylated protein
        p0 = 300;%[nM[ total concentration of protein (all unphosphorylated)
        V_kin_max = 10;%[nm/s] rate of activation (phosphorylation)
        V_phos_max = 10;%[nm/s] rate of de-phosphorylation
        KM1
        KM2
    end
    
    methods
        
       function self = MMBasicModelTwoProteinForms()
            self.KM1 = self.p0/2;
            self.KM2 = self.p0/2;
       end
       
       function [c,f,s] = pde_fun(self,x,t,u,DuDx)
            p = u(1);%unphosphorylated form of the protein
            pa = u(2);%active - phosphorylated form of the protein
            
            v_kin = self.V_kin_max * p ./ ( self.KM1 + p );%rate of activation
            v_kin = v_kin .* self.get_kinease_activity_normalizer(x,self.R);%effective rate of activation
            
            v_phos = self.V_phos_max * pa ./ ( self.KM2 + pa );%rate of deactivation
            
            c = [ 1 ; 1 ];
            f = [ self.D ; self.D ] .* DuDx;
            s = [ 2*self.D/x; 2*self.D/x ].*DuDx+[ -v_kin+v_phos ; +v_kin-v_phos ];
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

