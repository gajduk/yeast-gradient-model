classdef Cascade
    %CACSCADE - see Box 4 of Kholodenko signaling dyunamics in space and
    %time
    
    %Membrane          | 
    %                  V 
    %Cytoplasm    --| Mp1
    %                  | 
    %                  V 
    %             --| Mp2
    %                  | 
    %                  V 
    %             --| Mp3

    
    properties
        R = 10^-5;%cell radius [m]
        D = 5*10^-12;%diffusion constant
        
        %activation of M1 membrane
        v_mem_kin = 7.5;
        KM1 = 0.5;

        %DEactivation of M1 cytosolic
        V1_phos = 3;
        KMP1 = 0.7;
        
        %activation of M2 cytosolic
        k2_cat = 9;
        KM2 = 0.7;

        %DEactivation of M2 cytosolic
        V2_phos = 3;
        KMP2 = 0.7;
        %phosphorylation of M2 cytosolic
        k3_cat = 9;
        KM3 = 0.5;

        %DEactivation of M3 cytosolic
        V3_phos = 3;
        KMP3 = 0.7;
    end
        
    methods
        
       function self = Cascade()
       end
       
        function res = run(self,time_steps,end_time,spatial_steps)
            %parameters
            %time_steps,spatial_steps - number of discretization steps in time,space respectively (higher values better precision and
            %smoother solutions but longer running time
            %end_time - the system solutions at times 0 [s] - end_time [s] are
            %calculated
            xmesh = linspace(0,self.R,spatial_steps);
            tspan = linspace(0,end_time,time_steps);

            m = 2;
            
            pdefun = @(x,t,u,DuDx) self.pde_fun(x,t,u,DuDx);%actual PDE
            icfun = @(x) self.ic_fun(x);%initial conditions
            bcfun = @(xl,ul,xr,ur,t) self.bc_fun(xl,ul,xr,ur,t);%boundary conditions

            sol = pdepe(m,pdefun,icfun,bcfun,xmesh,tspan);
            
            res = kholodenko.KholodenkoPlotter(xmesh,tspan,sol);
        end
        
       function [v_mem_kin,v1_phos,v1_kin,v2_phos,v2_kin,v3_phos] = getRates(self,u)
            mp1 = u(1);
            mp2 = u(2);
            mp3 = u(3);
            m1 = 1-mp1;
            m2 = 1-mp2;
            m3 = 1-mp3;
            
            v_mem_kin = self.v_mem_kin .* m1 ./ (self.KM1+m1);
            v1_phos = self.V1_phos .* mp1 ./ (self.KMP1+mp1);
            
            v1_kin = self.k2_cat .* mp1 .* m2 ./ (self.KM2 + m2);
            v2_phos = self.V2_phos .* mp2 ./ (self.KMP2 + mp2);
            
            v2_kin = self.k3_cat .* mp2 .* m3 ./ (self.KM3 + m3);
            v3_phos = self.V3_phos .* mp3 ./ (self.KMP3 + mp3);
       end
       
       function [c,f,s] = pde_fun(self,x,t,u,DuDx)
            [~,v1_phos,v1_kin,v2_phos,v2_kin,v3_phos] = self.getRates(u);
            
            c = [1; 1; 1];
            f = [self.D; self.D; self.D] .* DuDx;
            s = [2*self.D/x; 2*self.D/x; 2*self.D/x] .* DuDx;
            
            s = s - [ +v1_phos ; 
                    +v2_phos - v1_kin ;
                    +v3_phos - v2_kin ];
       end
       
       function u0 = ic_fun(self,x)
            %u0 = [1; 1; 1];
            u0 = [0; 0; 0];
       end
       
       function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
            [v_mem_kin,~,~,~,~,~] = self.getRates(ur);
            % pl and ql correspond to the left boundary conditions (x = 0), and 
            % pr and qr correspond to the right boundary condition (x = R).

            pl = [0; 0; 0];
            ql = [1; 1; 1];
            pr = [ - v_mem_kin; 0; 0];
            qr = [1; 1; 1];
       end

    end
        
end

