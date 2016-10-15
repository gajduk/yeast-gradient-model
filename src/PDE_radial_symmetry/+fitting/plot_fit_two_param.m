function plot_fit_two_param(param1,param2)
close all
[WT,M] = utils.get_dataset();

[~,wtn] = size(WT);
[~,mn] = size(M);

WTmeans = mean(WT,2);
WTtemp = WT-min(WTmeans);
WTnorm = WTtemp./mean(WTtemp(1,:));
WTnorm_mean = mean(WTnorm,2);

Mmeans = mean(M,2);
Mtemp = M-min(WTmeans);
Mnorm = Mtemp./mean(Mtemp(1,:));
Mnorm_mean = mean(Mnorm,2);

[mr2, ~] = utils.rsquare(Mnorm,repmat(Mnorm_mean,1,mn));
[wtr2, ~] = utils.rsquare(WTnorm,repmat(WTnorm_mean,1,wtn));

fss=-7:.35:7;
n = length(fss);
res = cell(n,n);
r2wt = zeros(n,n);
r2m = zeros(n,n);
parfor i=1:n
    fs1 = fss(i);
    for k=1:n
        fs2 = fss(k);
        model = models.RealParam_YeastMM_ste11_ste7_fus3_ste5();
        param1_new_val = model.(param1)*2^fs1;
        model.(param1) = param1_new_val;
        param2_new_val = model.(param2)*2^fs2;
        model.(param2) = param2_new_val;
        WTf = utils.get_output(model);
        [r2wtA, ~] = utils.rsquare(WTnorm,repmat(WTf,1,wtn));
        [r2mA, ~] = utils.rsquare(Mnorm,repmat(WTf,1,mn));
        res{i,k} = WTf;
        r2wt(i,k) = r2wtA;
        r2m(i,k) = r2mA;
        i*100+k
    end
end
[FSSX,FSSY] = meshgrid(2.^fss,2.^fss);
figure('Position',[100,100,500,300])
contourf(FSSX,FSSY,r2wt./wtr2,'g');
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel(sprintf('Fold change in %s',param1));
ylabel(sprintf('Fold change in %s',param2));
xlim([2^-7,2^7])
ylim([2^-7,2^7])
c = colorbar;
title('WT')
c.Label.String = 'R^2 (normalized)';
export_fig(sprintf('WT_model_fits_two_param_%s_%s.png',param1,param2),'-transparent','-r600')


figure('Position',[100,100,500,300])
contourf(FSSX,FSSY,r2m./mr2,'g');
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel(sprintf('Fold change in %s',param1));
ylabel(sprintf('Fold change in %s',param2));
xlim([2^-7,2^7])
ylim([2^-7,2^7])
c = colorbar;
title('M')
c.Label.String = 'R^2 (normalized)';
export_fig(sprintf('M_model_fits_two_param_%s_%s.png',param1,param2),'-transparent','-r600')

end