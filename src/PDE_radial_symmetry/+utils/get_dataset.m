function [WT,M] = get_dataset()
WT = csvread('csv data filtered equally spaced\\Fus3 WT_filtered.txt');
M = csvread('csv data filtered equally spaced\\Fus3 Ste11-S243A_filtered.txt');

WT = WT(:,1:11)';
M = M(:,1:11)';
end

