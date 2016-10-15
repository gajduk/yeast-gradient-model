close all
figure;
[WT,M] = get_dataset();
x = 0:0.1:1;
Mmean = mean(M,2);
hold on
plot(x,Mmean,'r','Linewidth',3)
xlabel('Distance from shmoo tip')
ylabel('Alpha (Denoised)')
title('M')
plot(0:0.01:0.24,ones(1,25)*((Mmean(1)-Mmean(11))/2+Mmean(11)),'g')
xlim([0 1])
ylim([0 1])
grid on
plot(x,M,'b');
plot(ones(1,42)*0.24,0:0.01:0.41,'g')
legend({'Mean','Half gradient','Single cell'})
export_fig m_raw.png -transparent -r600

figure;
WTmean = mean(WT,2);
hold on
plot(x,WTmean,'r','Linewidth',3)
xlabel('Distance from shmoo tip')
ylabel('Alpha (Denoised)')
title('WT')
plot(0:0.01:0.29,ones(1,30)*((WTmean(1)-WTmean(11))/2+WTmean(11)),'g')
xlim([0 1])
ylim([0 1])
grid on
plot(x,WT,'b');
plot(ones(1,39)*0.285,0:0.01:0.38,'g')
legend({'Mean','Half gradient','Single cell'})
export_fig wt_raw.png -transparent -r600