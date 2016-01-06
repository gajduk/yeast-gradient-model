classdef ModelOutput < handle
    
    
    properties
        model = -1;
        xmesh = -1;
        tspan = -1;
        sol = -1;
        fraction_active_protein = -1;
        steady_state = -1;
        half_fraction = -1;
        R = -1;
    end
    
    methods
        function self = ModelOutput(model,xmesh,tspan,sol,fraction_active_protein,steady_state,half_fraction)
            self.model = model;
            self.xmesh = xmesh;
            self.tspan = tspan;
            self.sol = sol;
            self.fraction_active_protein = fraction_active_protein;
            self.steady_state = steady_state;
            self.half_fraction = half_fraction;
            self.R = model.R;
        end
        
        function plot_fraction(self)
            self.initfigure()
            surf(self.xmesh,self.tspan,self.fraction_active_protein,'EdgeColor','none');   
            self.correct_x_axis();
            model_meta = metaclass(self.model);
            title(model_meta.Name);
            xlabel('Distance from cell membrane [\mum]');
            ylabel('Time [s]');
            zlabel('Fraction of phosphorylated protein');
            colorbar;
        end
        
        function plot_fraction_steady_state(self)
            self.initfigure()
            plot(self.xmesh,self.fraction_active_protein(end,:));
            self.correct_x_axis();
            model_meta = metaclass(self.model);
            title(strcat(model_meta.Name,' - Steady state'));
            xlabel('Distance from cell membrane [\mum]');
            ylabel('Fraction of phosphorylated protein');
        end
        
        function plot_average_fraction(self)
            self.initfigure()
            average = squeeze(mean(self.sol,2));
            plot(self.tspan,average(:,2)./(average(:,2)+average(:,1)));
            model_meta = metaclass(self.model);
            title(strcat(model_meta.Name,' - Average fraction in the cell'));
            xlabel('Time [s]');
            ylabel('Fraction of phosphorylated protein');
        end
        
        function correct_x_axis(self)
            ticks = 0:self.R/10:self.R;
            tick_labels = cell(length(ticks),1);
            for i=1:length(ticks)
                tick_labels{length(ticks)-i+1} = sprintf('%.0f',ticks(i)/self.R*10);
            end
            set(gca,'Xtick',ticks,'XTickLabel',tick_labels)
        end
        
        function initfigure(self)
            figure('Position', [100, 100, 1500, 895]);
        end
    end
    
end

