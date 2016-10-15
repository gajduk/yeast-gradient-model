function plot_model_no_feedback()
close all
[WT,M] = utils.get_dataset();
plot_model_no_feedback_for(WT)
suptitle('WT')
export_fig wt_model_no_feedback.png -transparent -r600
plot_model_no_feedback_for(M)
suptitle('M')
export_fig m_model_no_feedback.png -transparent -r600

end

function plot_model_no_feedback_for(WT) 
[~,n] = size(WT);
x = 0:0.1:1;
WTmeans = mean(WT,2);
WTtemp = WT-min(WTmeans);
WTnorm = WTtemp./mean(WTtemp(1,:));
WTnorm_mean = mean(WTnorm,2);

figure('Position',[100 100 1000 500]);
subplot(1,2,1)
hold on
plot(x,WTnorm_mean,'r','LineWidth',2);
xlabel('Distance from shmoo')
ylabel('Activity (denoised and normalized)')
[r2m, ~] = utils.rsquare(WTnorm,repmat(WTnorm_mean,1,n));

model = models.RealParam_YeastMM_ste11_ste7_fus3_ste5();
model.ste11_act_vkin = 0;
WTf = utils.get_output(model);
plot(x,WTf,'g','LineWidth',2);
[r2f, ~] = utils.rsquare(WTnorm,repmat(WTf,1,n));
WTmean_hf = utils.get_hf(x,WTmeans);
WTf_hf = utils.get_hf(x,WTf);
legend({sprintf('Mean (R^2=%.3f, hf=%.3f)',r2m,WTmean_hf),sprintf('Model (R^2=%.3f, hf=%.3f)',r2f,WTf_hf)})
ylim([-.03,1.03])
subplot(1,2,2)
title(sprintf('Model (R^2=%.3f)',r2f))
utils.plot_residuals(x,WTf',WT,'g')
end
