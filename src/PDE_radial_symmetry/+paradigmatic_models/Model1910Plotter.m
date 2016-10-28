classdef Model1910Plotter
    properties
        titles = {'act kin freely diff','NON-act kin freely diff','act kin bound',...
                'NON-act kin bound','act. substr','NON-act substr'};
           
    end
    
    methods
        function plot_all_surf(self,xmesh,tspan,res)
             figure('Position',[100,0,1000,900]);
            k = 1;
            for i=[1,3,5,6]
                subplot(2,2,k)
                k = k+1;
                surf(xmesh,tspan,res(:,:,i),'EdgeColor','None')
                xlabel('x')
                set(gca,'XTick',[0,.0000099])
                set(gca,'XTickLabel',{'PM','NM'})
                ylabel('t')
                zlabel('Conc')
                title(self.titles{i});
            end
        end
        function plot_gradient(self,res,i,xmesh,label)
            y = res(end,xmesh>0.0000007,i);
            title_s = sprintf('%3.3f hf=%.3f',label,utils.get_hf(xmesh,y)/xmesh(end));
            plot(xmesh(xmesh>0.0000007),y,'DisplayName',title_s)
            set(gca,'XTick',[0.0000007,.0000099])
            xlim([0.0000007,.0000099])
            set(gca,'XTickLabel',{'PM','NM'})
            xlabel('x')
            ylabel('Conc')
            title(self.titles{i});
        end
    end
    
end

