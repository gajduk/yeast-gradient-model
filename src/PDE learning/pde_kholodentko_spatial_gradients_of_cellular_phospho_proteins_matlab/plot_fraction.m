function [sol,ka,kp] = plot_fraction(xmesh,tspan,sol,ka,kp)

u = sol(:,:,2)./(sol(:,:,2)+sol(:,:,1));
figure('Position', [100, 100, 1500, 895]);
surf(xmesh,tspan,u,'EdgeColor','none');    
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('Fraction of phosphorylated protein');
title(sprintf('kp=%0.1f   ka=%0.1f',kp,ka));
colorbar;

end

