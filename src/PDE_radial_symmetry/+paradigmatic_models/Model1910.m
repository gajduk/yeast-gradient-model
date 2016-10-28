classdef Model1910
    %Kinase has two forms, active and non-active; both of which bind to the
    %membrane with different rates
    %When the kinase is bound to the membrane is activates the substrate
    %The substrate has two forms, active and non_active. The active
    %substrate deactivates the kinase.

    properties
        kin_a0 = 100;%active kinase freely diffusing
        kin_n0 = 0;%NON-active kinase freely diffusing

        kin_a_b0 = 0;%active kinase bound to the membrane
        kin_n_b0 = 0;%NON-active kinase bound to the membrane

        subst_a0 = 0;%active freely diffusing substrate
        subst_n0 = 100;%NON-active freely diffusing substrate


        x_m = 0.05*10^-5;%this is the membrane

        D_subst = 1*10^-8;%rate of diffusion of the substrate
        
        D_kin = 5*10^-12;%rate of diffusion of the kinase

        b_on_a = .5;%binding on rate of active kinase
        b_on_n = .01;%binding on rate of NON-active kinase

        b_off_a = .5;%binding off rate of active kinase
        b_off_n = 1;%binding off rate of NON-active kinase

        s_Km = 50;%Km for activation of the substrate by the kinase
        s_k_cat = 1000;%k_cat rate of activation of the substrate by the kinase
        s_phos =  .1;%DE-activation rate of the substrate
        

        k_Km = 20;%Km for the DE-activation of the  by the substrate
        k_k_cat = 1;%k_cat rate of DE-activation of the kinase by the susbstrate
        k_act = 1;%activation rate of the kinase

        type = 1;
    end
    
    methods

        function obj = Model1910(type_)
            %type = 1, only kinase that is bound to the membrane can activate
            %the substrate
            %type = 2, only freely diffusing kinase can activate the substrate
            obj.type = type_;
        end
        
        function res = hill(self,u,K_m,coef)
            if nargin == 2
               coef = 1; 
            end
            res = u .^ coef ./ ( K_m + u.^coef );
        end

        function [c,f,s] = pde_fun(self,x,t,u,DuDx)
            kin_a = u(1);%active kinase freely diffusing
            kin_n = u(2);%NON-active kinase freely diffusing

            kin_a_b = u(3);%active kinase bound to the membrane
            kin_n_b = u(4);%NON-active kinase bound to the membrane

            subst_a = u(5);%active freely diffusing substrate
            subst_n = u(6);%NON-active freely diffusing substrate
            
            membr = x<self.x_m;
            b_on_a_xt = self.b_on_a .* membr .* kin_a;
            b_off_a_xt = self.b_off_a .* membr .* kin_a_b;
            
            b_on_n_xt = self.b_on_n .* membr .* kin_n;
            b_off_n_xt = self.b_off_n .* membr .* kin_n_b;

            if self.type == 1
                total_b_kin = kin_a_b+kin_n_b;
            else
                total_b_kin = kin_a;
            end
            
            
            s_v_kin = subst_n .* self.s_k_cat .* self.hill(total_b_kin,self.s_Km,1) .* 1000;
            s_v_phos = subst_a .* self.s_phos;
            
            k_v_act = kin_n .* self.k_act;
            k_v_deact = kin_a .* self.k_k_cat .* self.hill(subst_a,self.k_Km,1);


            b_on_n_xt = 0;
            b_off_n_xt = 0;
            k_v_act = 0;
            k_v_deact = 0;
            
            c = [1; 1; 1; 1; 1; 1];
            f = [self.D_kin; self.D_kin; 0; 0; self.D_subst; self.D_subst] .* DuDx;
            s = [2*self.D_kin/x; 2*self.D_kin/x; 0; 0; 2*self.D_subst/x; 2*self.D_subst/x;] .* DuDx;
            s = s+[-b_on_a_xt+b_off_a_xt+k_v_act-k_v_deact;
                   -b_on_n_xt+b_off_n_xt-k_v_act+k_v_deact;
                   +b_on_a_xt-b_off_a_xt;
                   +b_on_n_xt-b_off_n_xt;
                   +s_v_kin-s_v_phos;
                   -s_v_kin+s_v_phos; ];
        end

        function u0 = ic_fun(self,x)
            u0 = [self.kin_a0;self.kin_n0;
                self.kin_a_b0;self.kin_n_b0;
                self.subst_a0;self.subst_n0];
        end

        function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
            % pR is a constant and the initial concentration of phosphorylated protein at the membrane
            % D - diffusion constant
            % pl and ql correspond to the left boundary conditions (x = 0), and 
            % pr and qr correspond to the right boundary condition (x = R).
            pl = [0; 0; 0; 0; 0; 0;];
            ql = [1; 1; 1; 1; 1; 1;];
            pr = [0; 0; 0; 0; 0; 0;];
            qr = [1; 1; 1; 1; 1; 1;];
        end

        
        function [xmesh,tspan,res] = run(self,time_steps,end_time,spatial_steps)
            %parameters
            %time_steps,spatial_steps - number of discretization steps in time,space respectively (higher values better precision and
            %smoother solutions but longer running time
            %end_time - the system solutions at times 0 [s] - end_time [s] are
            %calculated
            R = 10^-5;%sphere radius
            xmesh = linspace(0,R,spatial_steps);
            tspan = linspace(0,end_time,time_steps);
        
            m = 2;
            pdefun = @(x,t,u,DuDx) self.pde_fun(x,t,u,DuDx);%actual PDE
            icfun = @(x) self.ic_fun(x);%initial conditions
            bcfun = @(xl,ul,xr,ur,t) self.bc_fun(xl,ul,xr,ur,t);%boundary conditions

            res = pdepe(m,pdefun,icfun,bcfun,xmesh,tspan);
        end
    end
end