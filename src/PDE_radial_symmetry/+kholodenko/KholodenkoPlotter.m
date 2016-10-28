classdef KholodenkoPlotter < handle
    
    properties
        titles = {'M1^p','M2^p','M3^p'};
        colors = {'green','blue', [1 0.5 0.2]};
        xmesh = -1;
        tspan = -1;
        sol = -1;
    end
    
    methods
        function self = KholodenkoPlotter(xmesh,tspan,sol)
            self.xmesh = xmesh;
            self.tspan = tspan;
            self.sol = sol;
        end
        
        function plot_all_surf(self)
            figure('Position',[100,100,1500,600]);
            k = 1;
            for i=[1,2,3]
                subplot(1,3,k)
                k = k+1;
                surf(self.xmesh,self.tspan,self.sol(:,:,i),'EdgeColor','None')
                xlabel('x')
                set(gca,'XTick',[0,.0000099])
                set(gca,'XTickLabel',{'PM','NM'})
                ylabel('t')
                zlim([0,1])
                zlabel('Conc')
                view(135, 30);
                title(self.titles{i});
            end
        end
        
        function plot_gradient(self,i,label)
            y = self.sol(end,:,i);
            title_s = sprintf('%s  [hf=%.3f]',label,utils.get_hf(self.xmesh,fliplr(y))/self.xmesh(end));
            plot(self.xmesh,fliplr(y),'DisplayName',title_s)
            set(gca,'XTick',linspace(0, self.xmesh(end),6))
            xlim([0,self.xmesh(end)])
            set(gca,'XTickLabel',{'PM','.2','.4','.6','.8','NM'})
            grid on
            xlabel('x')
            ylabel('Phos. fraction of kinase')
            ylim([0,1])
            title(self.titles{i});
        end
        
        function plot_gradients(self)
            figure('Position',[100,100,500,400]);
            for i=[1,2,3]
                y = self.sol(end,:,i);
                self.plot_gradient(i,self.titles{i})
                hold on
            end
            legend show
            title('Beautiful gradinets')
        end
    end
    
end

