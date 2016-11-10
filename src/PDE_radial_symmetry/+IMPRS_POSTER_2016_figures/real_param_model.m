function res = real_param_model()
    res = solve_pde();
end
function res = hill(x,Km,coef)
    %res = x.^coef ./ (x.^coef+Km);
    res = (tanh((x-Km)./coef)+1)/2;
end
function [c,f,s] = pde_fun(x,t,u,DuDx)
    R = 3*10^-6;%sphere radius [m] RealParam
    D = 4*10^-13;%diffusion constant of the phosphorylated protein RealParam

    ste11_0 = 40;%initial concentration [nM] RealParam
    fus3_0 = 200;%initial concentration [nM] RealParam
    ste5_0 = 40;%initial concentration [nM] RealParam

    alpha_length = 1000;%length of stimulation in time units

    alpha = 1;%binding on rate in precense of pheromone RealParam

    b_on = 0.011;%binding on rate of ste11 to ste5 RealParam
    b_off = 0.011;%binding off rate of ste11 to ste5

    b_on_p = 0.0001;%binding on rate of ste11_p to ste5
    b_off_p = 0.01;%binding off rate of ste11_p to ste5, should be equal to b_off


    fus3_act_vkin = 3.5;%activation of fus3 by ste5-ste11 complex
    fus3_act_KM = 15;%activation of fus3 by ste5-ste11 complex
    fus3_act_hill = 10;%hill coeficient
    phos_rate_fus3 = .7;%rate of de-phosphorylation of fus3 RealParam

    global ste11_act_vkin;
    %ste11_act_vkin = .1;%activation of ste11 by fus3p
    ste11_act_KM = 60;%activation of ste11 by fus3p
    ste11_act_hill = 3;%hill coeficient
    phos_rate_ste11 = .01;%rate of de-phosphorylation of ste11 RealParam
    
    if t > alpha_length
       alpha = 0; 
    end
            
    ste5 = u(1);%scaffold

    ste11 = u(2);%UNphosphorylated ste11
    ste11p = u(3);%phosphorylated ste11

    ste5_ste11 = u(4);%complex
    ste5_ste11p = u(5);%complex

    fus3 = u(6);%UNphosphorylated fus3
    fus3p = u(7);%phosphorylated fus3
    
    
    b_on_ = b_on * alpha * ste5 * ste11;
    b_off_ = b_off * ste5_0 * ste5_ste11; 

    b_on_p_ = b_on_p * alpha * ste5 * ste11p;
    b_off_p_ = b_off_p * ste5_0  * ste5_ste11p;
    
    ste5_cmplx = ste5_ste11+ste5_ste11p;
            
    fus3_act = fus3 * fus3_act_vkin * ...
        hill(ste5_cmplx,fus3_act_KM,fus3_act_hill);
    fus3_phos = phos_rate_fus3 * fus3p;

    ste11_act = ste11 * ste11_act_vkin * ...
        hill(fus3p,ste11_act_KM,ste11_act_hill);
    ste11_phos = phos_rate_ste11 *  ste11p;

    %ste11_act = 0;
    %ste11_phos = 0;
    
    c = [1; 1; 1; 1; 1; 1; 1];
    f = [0; D; D; 0; 0; D; D] .* DuDx;
    s = [0; 2*D/x; 2*D/x; 0; 0; 2*D/x; 2*D/x] .* DuDx;
    s = s+[ -b_on_-b_on_p_ + b_off_ + b_off_p_ ;%ste5
            -b_on_ + b_off_ - ste11_act + ste11_phos ;%ste11
            -b_on_p_ + b_off_p_ + ste11_act - ste11_phos ; %ste11p
            +b_on_-b_off_ ;% ste5_ste11
            +b_on_p_-b_off_p_ ;% ste5_ste11p
            -fus3_act + fus3_phos ;% fus3
            +fus3_act - fus3_phos];% fus3p

end
function u0 = ic_fun(x)
    R = 3*10^-6;%sphere radius [m] RealParam
    ste11_0 = 40;%initial concentration [nM] RealParam
    fus3_0 = 200;%initial concentration [nM] RealParam
    ste5_0 = 40;%initial concentration [nM] RealParam
        
    ste5 = ste5_0*(x>0.95*R);
    u0 = [ste5; ste11_0;  0;   0; 0; fus3_0; 0];
end
function [pl,ql,pr,qr] = bc_fun(xl,ul,xr,ur,t)
    % pR is a constant and the initial concentration of phosphorylated protein at the membrane
    % D - diffusion constant
    % pl and ql correspond to the left boundary conditions (x = 0), and 
    % pr and qr correspond to the right boundary condition (x = R).
    pl = [0; 0; 0; 0; 0; 0; 0];
    ql = [1; 1; 1; 1; 1; 1; 1];
    pr = [0; 0; 0; 0; 0; 0; 0];
    qr = [1; 1; 1; 1; 1; 1; 1];
end

function res = solve_pde()
    R = 3*10^-6;
    xmesh = linspace(0,R,100);
    tspan = linspace(0,60,100);
    m = 2;
    res = pdepe(m,@pde_fun,@ic_fun,@bc_fun,xmesh,tspan);
end