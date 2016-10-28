
classdef RealParam_YeastMM_ste11_ste7_fus3_ste5 < ModelCore
    
    %RealParam_YeastMM_ste11_ste7_fus3_ste5
    
    properties
        R = 3*10^-6;%sphere radius [m] RealParam
        D = 4*10^-12;%diffusion constant of the phosphorylated protein RealParam
        
        ste11_0 = 40;%initial concentration [nM] RealParam
        fus3_0 = 200;%initial concentration [nM] RealParam
        ste5_0 = 40;%initial concentration [nM] RealParam
        
        phos_rate_fus3 = 10;%rate of de-phosphorylation of fus3 RealParam
        phos_rate_ste11 = 5;%rate of de-phosphorylation of ste11 RealParam
        
        alpha_length = 300;%length of stimulation in time units
        
        
        b_on = .1;%binding on rate of ste11 to ste5 RealParam
        alpha = .9%binding on rate in precense of pheromone RealParam
        b_off = 5;%binding off rate of ste11 to ste5
        
        %parameters to be fitted
        b_on_p = .01;%binding on rate of ste11_p to ste5
        b_off_p = 5;%binding off rate of ste11_p to ste5, should be equal to b_off
        
        fus3_act_vkin = 200;%activation of fus3 by ste5-ste11 complex
        fus3_act_KM = 10;%activation of fus3 by ste5-ste11 complex
        
        %parameters to be fitted
        ste11_act_vkin = 500;%activation of ste11 by fus3p
        ste11_act_KM = 20;%activation of ste11 by fus3p
    end
        
    methods
        
       function self = RealParam_YeastMM_ste11_ste7_fus3_ste5()
       end
       
       function [c,f,s] = pde_fun(self,x,t,u,DuDx)
            alpha = self.alpha;
            if t > self.alpha_length
               alpha = 0; 
            end
            
            ste5 = u(1);%scaffold
            
            ste11 = u(2);%UNphosphorylated ste11
            ste11p = u(3);%phosphorylated ste11
            
            ste5_ste11 = u(4);%complex
            ste5_ste11p = u(5);%complex
            
            fus3 = u(6);%UNphosphorylated fus3
            fus3p = u(7);%phosphorylated fus3
            
            b_on_ = ( self.b_on + alpha ) * ste5 * ste11;
            b_off_ = self.b_off * ste5_ste11;
            
            b_on_p_ = self.b_on_p * ste5 * ste11p;
            b_off_p_ = self.b_off_p * ste5_ste11p;
            
            fus3_act = fus3 * self.fus3_act_vkin * ( (ste5_ste11+ste5_ste11p) / ( self.fus3_act_KM + (ste5_ste11+ste5_ste11p) ) );
            fus3_phos = self.phos_rate_fus3 * fus3p;
            
            ste11_act = ste11 * self.ste11_act_vkin * ( fus3p / (self.ste11_act_KM+fus3p) );
            ste11_phos = self.phos_rate_ste11 *  ste11p;
            
            c = [1; 1; 1; 1; 1; 1; 1];
            f = [0; self.D; self.D; 0; 0; self.D; self.D] .* DuDx;
            s = [0; 2*self.D/x; 2*self.D/x; 0; 0; 2*self.D/x; 2*self.D/x] .* DuDx;
            s = s+[ -b_on_-b_on_p_ + b_off_ + b_off_p_ ;%ste5
                    -b_on_ + b_off_ - ste11_act + ste11_phos ;%ste11
                    -b_on_p_ + b_off_p_ + ste11_act - ste11_phos ; %ste11p
                    +b_on_-b_off_ ;% ste5_ste11
                    +b_on_p_-b_off_p_ ;% ste5_ste11p
                    -fus3_act + fus3_phos ;% fus3
                    +fus3_act - fus3_phos];% fus3p
                    
       end
       
       function u0 = ic_fun(self,x)
           
            u0 = [self.ste5_0*(x>self.R*0.95); self.ste11_0; 0; 0; 0; self.fus3_0; 0];
       end
       
       function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
            % pR is a constant and the initial concentration of phosphorylated protein at the membrane
            % D - diffusion constant
            % pl and ql correspond to the left boundary conditions (x = 0), and 
            % pr and qr correspond to the right boundary condition (x = R).
            pl = [0; 0; 0; 0; 0; 0; 0];
            ql = [1; 1; 1; 1; 1; 1; 1];
            pr = [0; 0; 0; 0; 0; 0; 0];
            qr = [1; 1; 1; 1; 1; 1; 1];
       end

    end
    
end

