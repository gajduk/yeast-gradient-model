function plot_residuals(x,fit,data,color)
mu = mean(data,2)';
sigma = std(data');
y = (fit-mu)./sigma;
yh = zeros(size(x));
plot(x,y,strcat(color,'o'))
hold on
plot(x,yh,'r')
arrayfun(@(x,y,yh) line([x x],[y yh],'color',color),x,y,yh)
xlim([-0.05 1.05])
end

