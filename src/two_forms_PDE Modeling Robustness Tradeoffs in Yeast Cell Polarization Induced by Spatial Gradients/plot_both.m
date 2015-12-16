function sol = plot_both(xmesh,tspan,sol)
u1 = sol(:,:,1);
u2 = sol(:,:,2);
lim = [0 1];
figure('Position', [100, 100, 1500, 895]);
ax(1) = subplot(1,2,1);
surf(xmesh,tspan,u1)    
xlabel('Distance x [m]')
ylabel('Time t [s]')
zlabel('Concentration of UNphosphorylated protein')
ax(2) = subplot(1,2,2);
surf(xmesh,tspan,u2)    
xlabel('Distance x [m]')
ylabel('Time t [s]')
zlabel('Concentration of phosphorylated protein')
h=colorbar;
set(h, 'Position', [.9 .11 .015 .8150])
caxis(lim);
set(ax,'CLim',lim);
for i=1:2
      pos=get(ax(i), 'Position');
      set(ax(i), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
end
end

