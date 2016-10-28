function res = foling_around()
    solve_pde();
end

function [c,f,s] = pde_fun(x,t,u,DuDx)
    r1 = 1;
    r2 = 1;
    r3 = 1;
    
    h1 = 3;
    h2 = 3;
    h3 = 3;
    
    c = [1; 1; 1;];
    D = 10;
    f = [1; 1; 1; ] * D .* DuDx;
    
    s = [r1 * ( (1-u(1))-u(1) * 10*u(3).^h1./(1+u(3).^h1) ) ; ...
        r2 * ( (1-u(2)) * u(1).^h2./(1+u(1).^h2)-.1*u(2) ) ; ...
        r3 * ( (1-u(3)) * u(2).^h3./(1+u(2).^h3)-.1*u(3) ) ];
end
function u0 = ic_fun(x)
    u0 = [0;0;0];
            
end
function [pl,ql,pr,qr] = bc_fun(xl,ul,xr,ur,t)
    % pR is a constant and the initial concentration of phosphorylated protein at the membrane
    % D - diffusion constant
    % pl and ql correspond to the left boundary conditions (x = 0), and 
    % pr and qr correspond to the right boundary condition (x = R).
    pl = [0; 0; 0;];
    ql = [1; 1; 1;];
    pr = [0; 0; 0;];
    qr = [1; 1; 1;];
end

function res = solve_pde()
        xmesh = linspace(0,10,100);
    	tspan = linspace(0,100,100);
        m = 2;
        res = pdepe(m,@pde_fun,@ic_fun,@bc_fun,xmesh,tspan);
        figure('Position',[100,100,1000,300]);
        subplot(1,3,1)
        surf(xmesh,tspan,res(:,:,1),'EdgeColor','None')
        xlabel('x')
        ylabel('t')
        subplot(1,3,2)
        surf(xmesh,tspan,res(:,:,2),'EdgeColor','None')
        xlabel('x')
        ylabel('t')
        
        subplot(1,3,3)
        surf(xmesh,tspan,res(:,:,3),'EdgeColor','None')
        xlabel('x')
        ylabel('t')
end
function res = solve_ode()

    x0 = [0;0;0];
    res = ode45(@(t,x) ode_fun(t,x),[0 100],x0);
    plot(res.x,res.y);
    legend({'1','2','3'});
end

function dx = ode_fun(t,x)
    r1 = 1;
    r2 = .5;
    r3 = .5;
    
    h1 = 3;
    h2 = 3;
    h3 = 3;
    dx = [r1 * ( (1-x(1))-x(1) * 10*x(3).^h1./(1+x(3).^h1) ) ; ...
        r2 * ( (1-x(2)) * x(1).^h2./(1+x(1).^h2)-.1*x(2) ) ; ...
        r3 * ( (1-x(3)) * x(2).^h3./(1+x(2).^h3)-.1*x(3) ) ];
end
