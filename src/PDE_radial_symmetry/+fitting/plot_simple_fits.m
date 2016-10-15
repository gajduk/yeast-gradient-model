function plot_simple_fits()
close all
[WT,M] = utils.get_dataset();
figure('Position',[100 100 1000 500]);
plot_simple_fits_for(M)
export_fig m_simple_fit.png -transparent -r600
figure('Position',[100 100 1000 500]);
plot_simple_fits_for(WT)
export_fig wt_simple_fit.png -transparent -r600
end

function plot_simple_fits_for(M)
[~,n] = size(M);
subplot(2,4,[1,2,5,6]);
x = 0:0.1:1;
Mmean = mean(M,2);
Msigma = std(M');
plot(x,Mmean,'r','Linewidth',3)
[r2m ~] = utils.rsquare(M,repmat(Mmean,1,n));
hold on
xx = repmat(x,n,1);
p1 = polyfit(xx,M',1);
M1 = polyval(p1,x);
[r21 ~] = utils.rsquare(M',repmat(M1,n,1));
p2 = polyfit(xx,M',2);
M2 = polyval(p2,x);
[r22 ~] = utils.rsquare(M',repmat(M2,n,1));
p3 = polyfit(xx,M',3);
M3 = polyval(p3,x);
[r23 ~] = utils.rsquare(M',repmat(M3,n,1));
MM = M';
[f,gof,~] =  fit(xx(:),MM(:),'Exp1');
plot(x,f(x),'g-')
plot(x,M1,'c-','Linewidth',1)
plot(x,M2,'b-','Linewidth',1)
plot(x,M3,'k-.','Linewidth',1)
xlabel('Distance from shmoo tip')
ylabel('Alpha (Denoised)')
title('M')
grid on
legend({sprintf('Mean (R^2=%.3f)',r2m),sprintf('Exp (R^2=%.3f)',gof.rsquare),sprintf('Fit poly 1 (R^2=%.3f)',r21),sprintf('Fit poly 2 (R^2=%.3f)',r22),sprintf('Fit poly 3 (R^2=%.3f)',r23)})
subplot(2,4,3);
utils.plot_residuals(x,M1,M,'c');
title(sprintf('Fit poly 1 (R^2=%.3f)',r21))
subplot(2,4,4);
utils.plot_residuals(x,M2,M,'b');
title(sprintf('Fit poly 2 (R^2=%.3f)',r22))
subplot(2,4,7);
utils.plot_residuals(x,M3,M,'k');
title(sprintf('Fit poly 3 (R^2=%.3f)',r23))
subplot(2,4,8);
utils.plot_residuals(x,f(x)',M,'g');
title(sprintf('Exp (R^2=%.3f)',gof.rsquare))
end