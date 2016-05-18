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
        time_scale = 's';
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
            if self.tspan(end) > 1000
               self.tspan = tspan./60;
               self.time_scale = 'min';
            end
        end
        
        function plot_all(self)
           set(0,'DefaultAxesFontSize',9)
           set(0,'DefaultTextFontSize',9)
           figure('Position', [100, 100, 900, 196]);
           subplot(1,3,1)
           self.plot_fraction_internal()
           title('')
           subplot(1,3,2)
           self.plot_fraction_steady_state_internal()
           title('Steady state')
           subplot(1,3,3)
           self.plot_average_fraction_internal()
           title('Cell average')
           export_fig(self.get_model_name(),'-png','-transparent','-r300')
        end
        
        function plot_fraction(self)
           self.initfigure()
           self.plot_fraction_internal()
        end
        
        function model_name = get_model_name(self)
            model_meta = metaclass(self.model);
            model_name = model_meta.Name;
        end
        
        function plot_phosphatase_internal(self)
            surf(self.xmesh,self.tspan,fliplr(self.sol(:,:,4)./(self.sol(:,:,3)+self.sol(:,:,4))),'EdgeColor','none');   
            self.correct_x_axis();
            title(self.get_model_name());
            xlabel({'Distance from','plasma membrane [\mum]'});
            ylim([0 max(self.tspan)])
            ylabel(strcat('Time [',self.time_scale,']'));
            zlabel('Fraction of active phosphatase');
            %colorbar;
            view([50 52]);
        end
                
        function plot_all_sol(self)
            [~,~,sols_len] = size(self.sol);
            figure;
            for idx=1:sols_len
                subplot(1,sols_len,idx);
                surf(self.xmesh,self.tspan,fliplr(self.sol(:,:,idx)),'EdgeColor','none');   
                self.correct_x_axis();
                title(self.get_model_name());
                xlabel({'Distance from','plasma membrane [\mum]'});
                ylim([0 max(self.tspan)])
                ylabel(strcat('Time [',self.time_scale,']'));
                zlabel('Fraction of active protein');
            %colorbar;
                view([50 52]);
            end
        end
        
        
        function plot_fraction_internal(self)
            surf(self.xmesh,self.tspan,fliplr(self.fraction_active_protein),'EdgeColor','none');   
            self.correct_x_axis();
            title(self.get_model_name());
            xlabel({'Distance from','plasma membrane [\mum]'});
            ylim([0 max(self.tspan)])
            ylabel(strcat('Time [',self.time_scale,']'));
            zlabel('Fraction of active protein');
            %colorbar;
            view([50 52]);
        end
        
        function plot_fraction_steady_state(self)
            self.initfigure()
            self.plot_fraction_steady_state_internal()
        end
            
        function plot_fraction_steady_state_internal(self,normalized_param)
            normalized = false;
            if nargin > 1
                normalized = normalized_param;
            end
            y = fliplr(self.fraction_active_protein(end,:));
            if normalized
                y = y./self.fraction_active_protein(end,end);
            end
            plot(self.xmesh,y);
            self.correct_x_axis();
            title(strcat(self.get_model_name(),' - Steady state'));
            xlabel('Distance from plasma membrane [\mum]');
            if normalized
                ylabel({'Fraction of active protein','Normalized to fraction at PM'});
            else
                ylabel('Fraction of active protein');
            end
            grid on
        end
        
        function plot_average_fraction(self)
            self.initfigure()
            self.plot_average_fraction_internal()
        end
        
        function plot_average_fraction_internal(self,normalized_param)
            normalized = false;
            if nargin > 1
                normalized = normalized_param;
            end
            average = squeeze(mean(self.sol,2));
            y = average(:,2)./(average(:,2)+average(:,1));
            if normalized
                y = y./y(end);
            end
            plot(self.tspan,y);
            title(strcat(self.get_model_name(),' - Average fraction in the cell'));
            xlabel(strcat('Time [',self.time_scale,']'));
            ylabel('Fraction of active protein');
            if normalized
                ylabel({'Fraction of active protein','Normalized to avg at SS'});
            else
                ylabel('Fraction of active protein');
            end
            xlim([0 max(self.tspan)])
            grid on
        end
        
        function correct_x_axis(self)
            ticks = 0:self.R/5:self.R;
            tick_labels = cell(length(ticks),1);
            for i=1:length(ticks)
                tick_labels{i} = sprintf('%.0f',ticks(i)/self.R*10);
            end
            tick_labels{1} = 'PM';
            tick_labels{length(tick_labels)} = 'NM';
            set(gca,'XTick',ticks,'XTickLabel',tick_labels)
        end
        
        function res = highest_fraction(self)
            average = squeeze(mean(self.sol,2));
            y = average(:,2)./(average(:,2)+average(:,1));
            y = y./y(end);
            res = max(y);
        end
        
        function res = transient_time(self)
            average = squeeze(mean(self.sol,2));
            y = average(:,2)./(average(:,2)+average(:,1));
            res = self.tspan(y==max(y));
        end
                
        function initfigure(self)
            figure('Position', [100, 100, 1500, 895]);
            
        end
    end
    
end

