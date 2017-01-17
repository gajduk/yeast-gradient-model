
classdef RealParam_YeastMM_ste11_ste7_fus3_ste5 < ModelCore
    
    %RealParam_YeastMM_ste11_ste7_fus3_ste5
    
    properties
        R = 3*10^-6;%sphere radius [m] RealParam
        D = 4*10^-12;%diffusion constant of the phosphorylated protein RealParam
        
        ste11_0 = 40;%initial concentration [nM] RealParam
        fus3_0 = 200;%initial concentration [nM] RealParam
        ste5_0 = 40;%initial concentration [nM] RealParam
        
        alpha_length = 1000;%length of stimulation in time units
        
        alpha = 1;%binding on rate in precense of pheromone RealParam
        
        b_on = 0.1;%binding on rate of ste11 to ste5 RealParam
        b_off = 0.1;%binding off rate of ste11 to ste5
        
        b_on_p = 0.008;%binding on rate of ste11_p to ste5
        b_off_p = 0.1;%binding off rate of ste11_p to ste5, should be equal to b_off
        
        
        fus3_act_vkin = 10;%activation of fus3 by ste5-ste11 complex
        fus3_act_KM = 10;%activation of fus3 by ste5-ste11 complex
        fus3_act_hill = 10;%hill coeficient
        phos_rate_fus3 = 10;%rate of de-phosphorylation of fus3 RealParam
       
        
        ste11_act_vkin = 10;%activation of ste11 by fus3p
        ste11_act_KM = .2;%activation of ste11 by fus3p
        ste11_act_hill = 5;%hill coeficient
        
        phos_rate_ste11 = .3;%rate of de-phosphorylation of ste11 RealParam
       
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
            
            b_on_ = self.b_on * alpha * ste5 * ste11;
            b_off_ = self.b_off * self.ste5_0 * ste5_ste11; 
            
            b_on_p_ = self.b_on_p * alpha * ste5 * ste11p;
            b_off_p_ = self.b_off_p * self.ste5_0  * ste5_ste11p;
            
            ste5_cmplx = ste5_ste11+ste5_ste11p;
            
            fus3_act = fus3 * self.fus3_act_vkin * ...
                self.hill(ste5_cmplx/self.ste5_0,self.fus3_act_KM,self.fus3_act_hill);
            fus3_phos = self.phos_rate_fus3 * fus3p;
            
            ste11_act = ste11 * self.ste11_act_vkin * ...
                self.hill(fus3p/self.fus3_0,self.ste11_act_KM,self.ste11_act_hill);
            ste11_phos = self.phos_rate_ste11 *  ste11p;
            
            ste11_act = 0;
            ste11_phos = 0;
            
            c = [1; 1; 1; 1; 1; 1; 1; 1];
            f = [0; self.D; self.D; 0; 0; self.D; self.D; 0] .* DuDx;
            s = [0; 2*self.D/x; 2*self.D/x; 0; 0; 2*self.D/x; 2*self.D/x; 0] .* DuDx;
            s = s+[ -b_on_-b_on_p_ + b_off_ + b_off_p_ ;%ste5
                    -b_on_ + b_off_ - ste11_act + ste11_phos ;%ste11
                    -b_on_p_ + b_off_p_ + ste11_act - ste11_phos ; %ste11p
                    +b_on_-b_off_ ;% ste5_ste11
                    +b_on_p_-b_off_p_ ;% ste5_ste11p
                    -fus3_act + fus3_phos ;% fus3
                    +fus3_act - fus3_phos;% fus3p
                    self.extra_fun(u)-u(8) ];
                    
       end
       
       function res = extra_fun(self,u)
            ste5 = u(1);%scaffold
            
            ste11 = u(2);%UNphosphorylated ste11
            ste11p = u(3);%phosphorylated ste11

            ste5_ste11 = u(4);%complex
            ste5_ste11p = u(5);%complex

            fus3 = u(6);%UNphosphorylated fus3
            fus3p = u(7);%phosphorylated fus3
            
            ste5_cmplx = ste5_ste11+ste5_ste11p;
            
            res = ste5_cmplx/self.ste5_0; 
       end
       
       function res = hill(self,x,Km,coef)
           res = x.^coef ./ ( Km + x.^coef);
       end
       
       function u0 = ic_fun(self,x)
            ste5 = self.ste5_0*(x>0.95*self.R);
            u0 = [ste5; self.ste11_0;  0;   0; 0; self.fus3_0; 0];
            u0 = [u0;self.extra_fun(u0)];
       end
       
       function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
            % pR is a constant and the initial concentration of phosphorylated protein at the membrane
            % D - diffusion constant
            % pl and ql correspond to the left boundary conditions (x = 0), and 
            % pr and qr correspond to the right boundary condition (x = R).
            pl = [0; 0; 0; 0; 0; 0; 0; 0];
            ql = [1; 1; 1; 1; 1; 1; 1; 1];
            pr = [0; 0; 0; 0; 0; 0; 0; 0];
            qr = [1; 1; 1; 1; 1; 1; 1; 1];
       end

    end
    
end

