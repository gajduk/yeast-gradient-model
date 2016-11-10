classdef BindingPlotter < handle
    
    properties
        titles = {'k','k_b','s'};
        colors = {'green','blue', [1 0.5 0.2]};
        xmesh = -1;
        tspan = -1;
        sol = -1;
    end
    
    methods
        function self = BindingPlotter(xmesh,tspan,sol)
            self.xmesh = xmesh;
            self.tspan = tspan;
            self.sol = sol;
        end
        
        function plot_all_surf(self)
            figure('Position',[100,100,1500,600]);
            
            for i=[1,2,3]
                subplot(1,3,i)
                if i == 2
                    plot(self.tspan,self.sol(:,end,i))
                    xlabel('t')
                    ylabel('Conc')
                    ylim([0,1])
                else
                    surf(self.xmesh,self.tspan,self.sol(:,:,i),'EdgeColor','None')
                    zlim([0,1])
                    zlabel('Conc')
                    view(135, 30);
                    ylabel('t')
                    xlabel('x')
                    set(gca,'XTick',[0,.0000099])
                    set(gca,'XTickLabel',{'NM','PM'})
                end
                title(self.titles{i});
            end


        end
        
        function plot_gradient(self,i,label,norma)
            if nargin == 3
                norma = false;
            end

            y = self.sol(end,:,i);
            if norma
                y = y-min(y);
                y = y/max(y);
            end
            title_s = sprintf('%s  [hf=%.3f]',label,utils.get_hf(self.xmesh,fliplr(y))/self.xmesh(end));
            plot(self.xmesh,fliplr(y),'DisplayName',title_s)
            set(gca,'XTick',linspace(0, self.xmesh(end),6))
            xlim([0,self.xmesh(end)])
            set(gca,'XTickLabel',{'NM','.2','.4','.6','.8','PM'})
            grid on
            xlabel('x')
            ylabel('Phos. fraction of kinase')
            ylim([0,1])
            title(self.titles{i});
        end
        
        function plot_gradients(self)
            figure('Position',[100,100,700,400]);
            k = 1;
            for i=[1,3]
                subplot(1,2,k)
                k = k +1;
                y = self.sol(end,:,i);
                self.plot_gradient(i,self.titles{i})
                
                legend show
            end
            suptitle('Beautiful gradients')
        end
    end
    
end

