function fig3()
spatial_steps = 31;
time_steps = 100;
end_time = 30;
global ste11_act_vkin;
ste11_act_vkin = .1;
res1 = figures.real_param_model();

ste11_act_vkin = 0;
res2 = figures.real_param_model();

x = linspace(0,30,100);
figure('Position',[100 100 400 250])
plot(x,fliplr(res1(:,end,7)),'k',x,fliplr(res2(:,end,7)),'r')

xlabel('Time [min]')
ylabel('ppFus3')
legend({'WT','Ste11(S243A)'},'Location','south')
export_fig fig3.eps
end

