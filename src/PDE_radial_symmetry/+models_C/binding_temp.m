function binding_temp()
    feedback_analysis()
    export_fig('feedback_analysis.png')
end

function phos_analysis()
    global v_phos;
    global A;
    
    A = 2;
    
    v_phos0 = .01;
    
    figure('Position',[100,100,1000,300]);
    
    xmesh = linspace(0,1,100);
        
    for i=-4:2:4
        k = 2.^i;
        v_phos = v_phos0*k;
        res = solve_pde();
        
        if true
            
            subplot(1,2,1)
            norm(res(end,:,2))
            plot(xmesh,fliplr(res(end,:,2)),'DisplayName',sprintf('v_{phos} = %.2f * %.2f',v_phos0,k))
            xlabel('x')
            ylabel('Substrate')
            set(gca,'XTick',[0,1])
            set(gca,'XTickLabel',{'PM','NM'})
            hold on

            subplot(1,2,2)
            norm(res(end,:,2))
            plot(xmesh,normA(res(end,:,2)),'DisplayName',sprintf('v_{phos} = %.2f * %.2f',v_phos0,k))
            xlabel('x')
            ylabel('Substrate (norm)')
            set(gca,'XTick',[0,1])
            set(gca,'XTickLabel',{'PM','NM'})
            hold on
        end
    end
    subplot(1,2,1);
    legend('show');
    subplot(1,2,2);
    legend('show');
end

function feedback_analysis()
    global A;
    
    global v_phos;
    v_phos = 0.01;
    A0 = 16;
    
    figure('Position',[100,100,1000,300]);
    
    xmesh = linspace(0,1,100);
        
    for i=-4:2:4
        k = 2.^i;
        A = A0*k;
        res = solve_pde();
        
        if true
            
            subplot(1,2,1)
            norm(res(end,:,2))
            plot(xmesh,fliplr(res(end,:,2)),'DisplayName',sprintf('A = %.2f * %.2f',A0,k))
            xlabel('x')
            ylabel('Substrate')
            set(gca,'XTick',[0,1])
            set(gca,'XTickLabel',{'PM','NM'})
            hold on

            subplot(1,2,2)
            norm(res(end,:,2))
            plot(xmesh,normA(res(end,:,2)),'DisplayName',sprintf('A = %.2f * %.2f',A0,k))
            xlabel('x')
            ylabel('Substrate (norm)')
            set(gca,'XTick',[0,1])
            set(gca,'XTickLabel',{'PM','NM'})
            hold on
        end
    end
    subplot(1,2,1)
    legend('show')
    subplot(1,2,2)
    legend('show')
end

function res = normA(x)
x = x-x(1);
x = x/x(end);
res = fliplr(x);
end


function [c,f,s] = pde_fun(x,t,u,DuDx)
    D = .001;
    
    sbst = u(2);
    
    global v_phos;
    
    c = [1; 1; 1 ];
    
    f = [1; 1; 0 ] * D .* DuDx;
    
    s = [0 ; -sbst*v_phos ; 0];
end
function u0 = ic_fun(x)
    u0 = [1;0;0];
            
end


function [pl,ql,pr,qr] = bc_fun(xl,ul,xr,ur,t)
    % pl and ql correspond to the left boundary conditions (x = 0), and 
    % pr and qr correspond to the right boundary condition (x = R).
    k = ur(1);
    sbst = ur(2);
    k_b = ur(3);
    global A;
    
    b_on = .0001;
    b_off = .0001;
    
    v_kin = .01;
    
    pl = [0; 0; 0;];
    ql = [1; 1; 0;];
    pr = [+ k * b_on * (1-k_b) * (1+sbst)/(1+A * sbst) - k_b * b_off ;
        - v_kin * k_b * (1-sbst) ;
        + k_b * b_off  - k * b_on * (1-k_b)  * (1+sbst)/(1+A * sbst)];
    qr = [1; 1; 1;];
end

function res = solve_pde()
        xmesh = linspace(0,1,100);
    	tspan = linspace(0,1000,200);
        m = 2;
        res = pdepe(m,@pde_fun,@ic_fun,@bc_fun,xmesh,tspan);
        if false
            figure('Position',[100,100,1000,300]);
            subplot(1,3,1)
            surf(xmesh,tspan,res(:,:,1),'EdgeColor','None')
            xlabel('x')
            ylabel('t')
            zlabel('Kinase cytosolic')

            subplot(1,3,2)
            surf(xmesh,tspan,res(:,:,2),'EdgeColor','None')
            xlabel('x')
            ylabel('t')
            zlabel('Substrate cytosolic')

            subplot(1,3,3)
            plot(tspan,res(:,end,3))
            xlabel('t')
            ylabel('Kinase bound')
        end
end