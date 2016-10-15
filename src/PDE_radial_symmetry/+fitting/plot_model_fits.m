function plot_model_fits()
close all
[WT,M] = get_dataset();
x = 0:0.1:1;
Mmeans = mean(M,2);
WTmeans = mean(WT,2);
Mnorm = M/Mmeans(1);
WTnorm = WT/WTmeans(1);


figure;
hold on
plot(x,WTmeans/WTmeans(1),'r','LineWidth',3);
xlabel('Distance from shmoo')
ylabel('Alpha (denoised)')
title('WT')
WTnorm_mean = mean(WTnorm,2);
[r2m rmsem] = rsquare(WTnorm,repmat(WTnorm_mean,1,18));

spatial_steps = 41;
time_steps = 100;

end_time = 30;
model = RealParam_YeastMM_ste11_ste7_fus3_ste5();

minr2 = 0.65;
min_ste11_act_vkin = -1;
min_WTf = -1;
ste11_act_vkin = model.ste11_act_vkin;
for feedback_strength=-22:24:10
    model.ste11_act_vkin = ste11_act_vkin*2^feedback_strength;
    output1 = model.run(time_steps,end_time,spatial_steps);
    WTf = fliplr(output1.fraction_active_protein(end,:));
    WTf = WTf./output1.fraction_active_protein(end,end);
    WTf = WTf(:,1:4:41);
    [r2 rmse] = rsquare(WTnorm',repmat(WTf,18,1));
    feedback_strength
    r2
    if minr2 < r2
        minr2 = r2;
        min_ste11_act_vkin = ste11_act_vkin*2^feedback_strength;
        min_WTf = WTf;
    end
end
plot(0:0.1:1,min_WTf,'g-.','LineWidth',3);

plot(x,WTnorm,'b');
legend({sprintf('Mean (R^2=%.3f)',r2m),sprintf('Model fit (R^2=%.3f)',minr2),'Single cells'})
export_fig wt_fit_model.png -transparent -r600

figure;
hold on
plot(x,Mmeans/Mmeans(1),'r','LineWidth',3);
xlabel('Distance from shmoo')
ylabel('Alpha (denoised)')
title('M')
Mnorm_mean = mean(Mnorm,2);
[r2m rmsem] = rsquare(Mnorm,repmat(Mnorm_mean,1,25));

spatial_steps = 41;
time_steps = 100;

end_time = 30;
model = RealParam_YeastMM_ste11_ste7_fus3_ste5();

minr2 = 0.50;
min_ste11_act_vkin = -1;
min_Mf = -1;
ste11_act_vkin = model.ste11_act_vkin;
for feedback_strength=-22:24:10
    model.ste11_act_vkin = ste11_act_vkin*2^feedback_strength;
    output1 = model.run(time_steps,end_time,spatial_steps);
    Mf = fliplr(output1.fraction_active_protein(end,:));
    Mf = Mf./output1.fraction_active_protein(end,end);
    Mf = Mf(:,1:4:41);
    [r2 rmse] = rsquare(Mnorm',repmat(Mf,25,1));
    feedback_strength
    r2
    if minr2 < r2
        minr2 = r2;
        min_ste11_act_vkin = ste11_act_vkin*2^feedback_strength;
        min_Mf = Mf;
    end
end
plot(0:0.1:1,min_Mf,'g-.','LineWidth',3);

plot(x,Mnorm,'b');
legend({sprintf('Mean (R^2=%.3f)',r2m),sprintf('Model fit (R^2=%.3f)',minr2),'Single cells'})
export_fig m_fit_model.png -transparent -r600


end

