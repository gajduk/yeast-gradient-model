classdef (Abstract) ModelCore
    %MODELCORE A superclass for all models. Defines 3 abstract methods that
    %must be implement by each model.
    
    properties
    end
    
    methods (Abstract)
       pde_fun(self,x,t,u,DuDx)
       ic_fun(self,x)
       bc_fun(self,xl,ul,xr,ur,t)
    end
    
    methods
        
        function res = get_kinease_activity_normalizer(self,x,R)
            %parameters
            %x - the distance from the center of the cell, you get this as an
            %argument in the pde_fun
            %R - the radious of the cell
            a = 50;
            e_ar = exp(a*x/R);
            lambda = a/(exp(a)-1);
            res = e_ar*lambda/2;
        end
        
        function res = run(self,time_steps,end_time,spatial_steps)
            %parameters
            %time_steps,spatial_steps - number of discretization steps in time,space respectively (higher values better precision and
            %smoother solutions but longer running time
            %end_time - the system solutions at times 0 [s] - end_time [s] are
            %calculated
            R = self.R;
            xmesh = linspace(0,R,spatial_steps);
            tspan = linspace(0,end_time,time_steps);

            m = self.get_m();
            pdefun = @(x,t,u,DuDx) self.pde_fun(x,t,u,DuDx);%actual PDE
            icfun = @(x) self.ic_fun(x);%initial conditions
            bcfun = @(xl,ul,xr,ur,t) self.bc_fun(xl,ul,xr,ur,t);%boundary conditions

            sol = pdepe(m,pdefun,icfun,bcfun,xmesh,tspan);
            if size(sol,3) == 7
               u =  sol(:,:,7)./(sol(:,:,6)+sol(:,:,7));
            else
                if size(sol,3) > 1
                    u = sol(:,:,2)./(sol(:,:,2)+sol(:,:,1));
                else
                    u = sol(:,:,1);
                end 
            end
            
            steady_state = u(end,:);
            half = 1;
            for k=1:length(xmesh)
                if steady_state(k)-steady_state(1) < (steady_state(end)-steady_state(1))/2
                    half = k;
                end
            end
            half_fraction = xmesh(half)/R;
            r = xmesh;
            t = tspan;
            res = ModelOutput(self,r,t,sol,u,steady_state,half_fraction);
        end
        
        function m = get_m(self)
           m = 2;
        end
    end
    
end

