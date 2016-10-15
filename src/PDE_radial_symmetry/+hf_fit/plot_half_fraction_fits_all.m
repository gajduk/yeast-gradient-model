function plot_half_fraction_fits_all()
params = utils.get_all_model_params();
for i=1:length(params)
    params{i}
    plot_half_fraction_fits_for(params{i});
end
end

function plot_half_fraction_fits_for(param)
close all
[WT,M] = utils.get_dataset();
x = 0:0.1:1;
Mmeans = mean(M,2);
WTmeans = mean(WT,2);
Mhf = utils.get_hf(x,Mmeans);
WThf = utils.get_hf(x,WTmeans);
fss=-7:.25:7;
n = length(fss);
res = cell(1,n);
hfs = zeros(1,n);
parfor i=1:n
    fs = fss(i);
    model = models.RealParam_YeastMM_ste11_ste7_fus3_ste5();
    param_new_val = model.(param)*2^fs;
    model.(param) = param_new_val;
    WTf = utils.get_output(model);
    hf(i) = utils.get_hf(x,WTf);
end
figure('Position',[100,100,400,300])
semilogx(2.^fss,hf,'b')
hold on
semilogx(2.^fss,ones(size(fss))*Mhf,'k')
semilogx(2.^fss,ones(size(fss))*WThf,'r')
legend({'model','WT','M'})
xlabel(sprintf('Fold change in %s',param));
ylabel('Diff. in half fraction')
xlim([2^-7,2^7])
ylim([0.2,0.32])
title(sprintf('hf fit %s',param))
export_fig(sprintf('hf_fit_%s.png',param),'-transparent','-r600')
end

