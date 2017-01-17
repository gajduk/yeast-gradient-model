classdef D2Model

    % Three variables,
    % kinase - k ; bound_kinase - k_b ; substrate - s
    % the kinase binds and unbinds to/from the membrane
    % kinase and substrate are freely diffusing
    % ONLY the bound kinase can activate the substrate
    
    %Membrane            | <------
    %                    V       |           
    %Cytoplasm          <--      |          
    %				k_b	---> k   |    
    %   
    %                |       |   |          
    %                |       |   |          
    %                ---------   |        
    %                    |       |  
    %                    V       |      
    %                    s  -----| 
    %      


    properties
        R = 10^-5;%cell radius [m]
        D = 5*10^-12;%diffusion constant
        
        %binding of k on/off membrane
        input = 2;
        b_on = 0.000001;
        b_off = 0.000001;

        %DEactivation of k cytosolic
        k_phos = .3;
        KMPk = 0.1;

        %activation of substrate cytosolic
        ks_cat = .09;
        KMs = 0.7;

        %DEactivation of substrate cytosolic
        s_phos = .3;
        KMPs = 0.7;

        %feedback from susbstrate to kinase binding
        A = 1;
        kd = .5;
    end   

    methods

    	function [v_b_on,v_b_off,v_phos_k,v_kin_s,v_phos_s] = getRates(self,u)
			k = u(1);
			k_b = u(2);
			s = u(3);

			v_b_on = self.b_on .* ( 1-k_b) .* self.input .* ( 1+ self.A * s / self.kd ) ./ ( 1+  s / self.kd);
			v_b_off = self.b_off .* k_b;

			v_phos_k = self.k_phos .* k ./ ( self.KMPk + k );

			v_phos_s = self.s_phos .* s ./ ( self.KMPs + s );

			total_k = k_b;
			
			v_kin_s = self.ks_cat .* (1-s) .* total_k ./ ( self.KMs + total_k);
		end
		function u0 = ic_fun(self,x)
		    u0 = [1;1;0];
		end   

		function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
		    % pl and ql correspond to the left boundary conditions (x = 0), and 
		    % pr and qr correspond to the right boundary condition (x = R).
		    [v_b_on,v_b_off,v_phos_k,v_kin_s,v_phos_s] = self.getRates(ur);
		    		    
		    pl = [0; 0; 0;];
		    ql = [1; 1; 1;];

		    pr = [ -v_b_off + v_b_on * ur(1) ; +v_b_off - v_b_on ; -v_kin_s];
		    qr = [1; 1; 1;];
		end  

		function [c,f,s] = pde_fun(self,x,t,u,DuDx)
		    [v_b_on,v_b_off,v_phos_k,v_kin_s,v_phos_s] = self.getRates(u);

		    c = [1; 1; 1 ];
		    
		    f = [1; 0; 1 ] * self.D .* DuDx;
		    
		    s = [-v_phos_k ; 0 ; v_kin_s-v_phos_s];
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
            
            res = models_D.DPlotter(xmesh,tspan,sol);
        end

	end    

end

