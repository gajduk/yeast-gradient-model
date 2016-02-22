
classdef AutocatalyticNegativeFeedbackTwoProteinForms < ModelCore
    %DirectNegativeFeedbackTwoProteinForms
    %See https://www.overleaf.com/2061351vhqrxw#/5221863/ - label phosphatase negative feedback
    properties (Constant)
        R = 1;%sphere radius
        D = 5*10^-3;%diffusion constant of the phosphorylated protein
        p0 = 1;%[nM]initial concentration of protein (all unphosphorylated)
        phos0 = 1;%[nM] initial concentration of phospatases (all unactive)
    end
    
    properties
        ka = 0;%s^-1 rate of activation (phosphorylation)
        kp = .5;%s^-1 rate of de-phosphorylation
        
        kac = 1.5;%s^-1 autocatalytic activation
        kan = 0;%s^-1 autonomous activation
        
        phos_ka = .3;%phosphataset activation
        phos_kd = .1;%phosphatase de-activation
        
    end
    
    methods
       function [c,f,s] = pde_fun(self,x,t,u,DuDx)
            p = u(1);%unphosphorylated form of the protein
            pa = u(2);%active - phosphorylated form of the protein
            phos = u(3);%non-active phosphatase
            phosa = u(4);%active phosphatase
            ea = self.ka*p;% * x^20 * (sin(t/3)+1);
            d = .2;%rate fo diffusion of protein compared to phosphatase
            
            c = [1; 1; 1; 1];
            f = [d*self.D; d*self.D; self.D; self.D].*DuDx;
            s = [d*2*self.D/x; d*2*self.D/x; 2*self.D/x; 2*self.D/x].*DuDx;
            %%     autonomous       autocatalytic    phosphorylation   activation
            s = s+[ -self.kan * p - self.kac * pa * p + self.kp*pa*phosa - ea;%non-active protein
                    +self.kan * p + self.kac * pa * p - self.kp*pa*phosa + ea;%active protein
              ... %  phos activation          phos DE-activation
                    -self.phos_ka*pa*phos + self.phos_kd*phosa;%phosphotases
                     self.phos_ka*pa*phos - self.phos_kd*phosa];%active phosphotases
       end
       
       function u0 = ic_fun(self,x)
           if x > self.R*0.48 && x < self.R*0.52
               p = 0.5+rand()/20;
               phos = 0.25+rand()/20;
               u0 = [self.p0*(1-p); self.p0*p; self.phos0*(1-phos); self.phos0*phos];
           else
               u0 = [self.p0; 0; self.phos0; 0];
           end
           %u0 = [self.p0; 0; self.phos0; 0];
       end
       
       function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
            % pl and ql correspond to the left boundary conditions (x = 0), and 
            % pr and qr correspond to the right boundary condition (x = R).
            pl = [0; 0; 0; 0];
            ql = [1; 1; 1; 1];
            pr = [0; 0; 0; 0];
            qr = [1; 1; 1; 1];
       end

    end
    methods(Static)
        function sensitivity_analysis_kp()
            close all
            time_steps = 2000;
            end_time = 7200;
            spatial_steps = 500;
            kps = 10.^(-4:.05:4);
            half_fractions = zeros(length(kps),1);
            highest_fraction = zeros(length(kps),1);
            transient_time = zeros(length(kps),1);
            temp = InDirectNegativeFeedbackTwoProteinForms();
            initial_kp = temp.kp;
            
            parfor i=1:length(kps)
                model = InDirectNegativeFeedbackTwoProteinForms();
                model.kp = initial_kp*kps(i);
                output = model.run(time_steps,end_time,spatial_steps);
                half_fractions(i) = output.half_fraction;
                highest_fraction(i) = output.highest_fraction();
                transient_time(i) = output.transient_time();
            end
            
            figure('Position', [100, 100, 800, 270]);
            subplot(1,3,1);
            semilogx(kps,half_fractions);
            xlabel({'Regulation rate strength','relative to default'})
            ylabel({'Half fraction at SS','higher values \Rightarrow steeper gradient'})
            xlim([min(kps),max(kps)])
            ax = gca;
            ax.XTick = [10^-4,10^-2,1,10^2,10^4];
            ax.XTickLabel = {'10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}'};
            grid on
            
            subplot(1,3,2);
            loglog(kps,highest_fraction);
            xlabel({'Regulation rate strength','relative to default'})
            ylabel({'Peak fraction','Normalized to avg at SS'})
            xlim([min(kps),max(kps)])
            ax = gca;
            ax.XTick = [10^-4,10^-2,1,10^2,10^4];
            ax.XTickLabel = {'10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}'};
            grid on
            
            subplot(1,3,3);
            loglog(kps,transient_time);
            xlabel({'Regulation rate strength','relative to default'})
            ylabel('Peak time [min]')
            xlim([min(kps),max(kps)])
            ax = gca;
            ax.XTick = [10^-4,10^-2,1,10^2,10^4];
            ax.XTickLabel = {'10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}'};
            grid on
            
            export_fig('sensitivity_analysis_kp1','-png','-transparent','-r300')
        end
        
        function sensitivity_analysis_kp_sample()
            close all
            time_steps = 1000;
            end_time = 7200;
            spatial_steps = 50;
            kps = [0.01,0.1,0.3,1,3,10,100];
            
            %run models for each value in kps
            res = cell(length(kps),1);
            temp = InDirectNegativeFeedbackTwoProteinForms();
            initial_kp = temp.kp;
            
            parfor i=1:length(kps)
                model = InDirectNegativeFeedbackTwoProteinForms();
                model.kp = initial_kp*kps(i);
                output = model.run(time_steps,end_time,spatial_steps);
                res{i} = output;
            end
            
            figure('Position', [100, 100, 650, 220]);
            %plot the gradients at steady state 
            subplot('Position',[0.1, 0.2,0.25,0.7]);
            legend_labels = cell(1,length(res));
            for i=1:length(res)
                res{i}.plot_fraction_steady_state_internal(true);
                legend_labels{i} = sprintf('k_d = %.2f k^0_d',kps(i));
                hold on
            end
            title('Steady state');
            
            %plot the time profiles
            subplot('Position',[0.45, 0.2,0.50,0.7]);
            for i=1:length(res)
                res{i}.plot_average_fraction_internal(true);
                hold on
            end
            xlim([0 30])
            title('Cell average');
            grid on
            
            legend(legend_labels,'Location','eastoutside');
            export_fig('sensitivity_analysis_kp_sample','-png','-transparent','-r300')
        end
    end
    
end

