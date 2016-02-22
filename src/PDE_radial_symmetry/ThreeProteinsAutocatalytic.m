
classdef ThreeProteinsAutocatalytic < ModelCore
    %DirectNegativeFeedbackTwoProteinForms
    %See https://www.overleaf.com/2061351vhqrxw#/5221863/ - label phosphatase negative feedback
    properties (Constant)
        R = 1;%sphere radius
        D = 5*10^-3;%diffusion constant of the phosphorylated protein
        p0 = 100;%[nM]initial concentration of protein (all unphosphorylated)
        phos0 = 300;%[nM] initial concentration of phospatases (all unactive)
    end
    
    properties
        ka = 2;%s^-1 rate of activation due to ligand (phosphorylation)
        kp = .25;%s^-1 rate of de-phosphorylation
        
        kac = 0;%s^-1 autocatalytic activation
        kan = 0;%s^-1 autonomous activation
        
        phos_ka = .025;%phosphataset activation
        phos_kd = .75;%phosphatase de-activation
        
    end
    
    methods
       function [c,f,s] = pde_fun(self,x,t,u,DuDx)
            raf = u(1);%unphosphorylated form of the protein
            praf = u(2);%active - phosphorylated form of the protein
            
            mek = u(3);
            pmek = u(4);
            
            erk = u(5);%non-active phosphatase
            perk = u(6);%active phosphatase
            d = .2;%rate fo diffusion of protein compared to phosphatase
            
            ea = 0;
            if t < 1
                ea = self.ka * raf;
            end
            
            c = [1; 1; 1; 1; 1; 1];
            f = [d*self.D; d*self.D; self.D; self.D; self.D; self.D].*DuDx;
            s = [d*2*self.D/x; d*2*self.D/x; 2*self.D/x; 2*self.D/x; 2*self.D/x; 2*self.D/x].*DuDx;
            %% ref       autonomous           autocatalytic            phosphorylation           activation
            s = s+[  -self.kan * raf    -  self.kac * praf * raf  +  self.kp * praf * perk  -  ea;%non-active protein
                     +self.kan * raf    +  self.kac * praf * raf  -  self.kp * praf * perk  +  ea;%active protein
          ... %% mek   activation                DE-activation
                    -self.phos_ka*mek*praf + self.phos_kd*pmek;%phosphotases
                     self.phos_ka*mek*praf - self.phos_kd*pmek;
          ... %% erk   activation                DE-activation      
                    -self.phos_ka*erk*pmek + self.phos_kd*perk;%phosphotases
                     self.phos_ka*erk*pmek - self.phos_kd*perk];%active phosphotases
                 
       end
       
       function u0 = ic_fun(self,x)
           if x > self.R*0.48 && x < self.R*0.52
               p = 0.5+rand()/20;
               phos = 0.25+rand()/20;
               u0 = [self.p0*(1-p); self.p0*p; self.phos0*(1-phos); self.phos0*phos; self.phos0*(1-phos); self.phos0*phos];
           else
               u0 = [self.p0; 0; self.phos0; 0; self.phos0; 0];
           end
            u0 = [self.p0; 0; self.phos0; 0; self.phos0; 0];
       end
       
       function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
            % pl and ql correspond to the left boundary conditions (x = 0), and 
            % pr and qr correspond to the right boundary condition (x = R).
            pl = [0; 0; 0; 0; 0; 0];
            ql = [1; 1; 1; 1; 1; 1];
            pr = [0; 0; 0; 0; 0; 0];
            qr = [1; 1; 1; 1; 1; 1];
       end

    end
    
end

