function plot_r2_fit_all_params()
    params = utils.get_all_model_params();
    for i=1:length(params)
        params{i}
        plot_r2_fit_single_param(params{i});
    end
end

function res = plot_r2_fit_single_param(param)
close all
[WT,M] = utils.get_dataset();

[~,wtn] = size(WT);
[~,mn] = size(M);

WTmeans = mean(WT,2);
WTtemp = WT-min(WTmeans);
WTnorm = WTtemp./mean(WTtemp(1,:));
WTnorm_mean = mean(WTnorm,2);

Mmeans = mean(M,2);
Mtemp = M-min(Mmeans);
Mnorm = Mtemp./mean(Mtemp(1,:));
Mnorm_mean = mean(Mnorm,2);

[mr2, ~] = utils.rsquare(Mnorm,repmat(Mnorm_mean,1,mn));
[wtr2, ~] = utils.rsquare(WTnorm,repmat(WTnorm_mean,1,wtn));

fss=-7:.25:7;
n = length(fss);
res = cell(1,n);
parfor i=1:n
    fs = fss(i);
    
    model = models.RealParam_YeastMM_ste11_ste7_fus3_ste5();
    param_new_val = model.(param)*2^fs;
    model.(param) = param_new_val;
    WTf = utils.get_output(model);
    [r2wt, ~] = utils.rsquare(WTnorm,repmat(WTf,1,wtn));
    [r2m, ~] = utils.rsquare(Mnorm,repmat(WTf,1,mn));
    res{i}.fit = WTf;
    res{i}.r2wt = r2wt;
    res{i}.r2m = r2m;
    res{i}.fs = fs;
end
r2wt = zeros(size(fss));
r2m = zeros(size(fss));
for i=1:n
    r2wt(i) = res{i}.r2wt;
    r2m(i) = res{i}.r2m;
end
figure('Position',[100,100,500,300])
semilogx(2.^fss,r2wt./wtr2,'g');
hold on
semilogx(2.^fss,r2m./mr2,'c');
legend({sprintf('WT (R^2=%.3f)',wtr2),sprintf('M (R^2=%.3f)',mr2)},'Location','eastoutside')
xlabel(sprintf('Fold change in %s',param));
ylabel('R^2 (normalized )')
xlim([2^-7,2^7])
ylim([0.92,1])
grid on
export_fig(sprintf('R2_fits_%s.png',param),'-transparent','-r600')
end

