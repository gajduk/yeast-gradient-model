
classdef PatternFormation < ModelCore
    %DirectNegativeFeedbackTwoProteinForms
    %See https://www.overleaf.com/2061351vhqrxw#/5221863/ - label phosphatase negative feedback
    properties (Constant)
        R = 10^-5;%sphere radius
        D = 5*10^-14;%diffusion constant of the phosphorylated protein
        p0 = 300;%[nM]initial concentration of protein (all unphosphorylated)
        phos0 = 50;%[nM] initial concentration of phospatases (all unactive)
        phos_feedback = .002;% - the speed of phosphatase dynamics
        A = 1;% strength of the feedback (a free parameter).
    end
    
    properties
    end
    
    methods
       function [c,f,s] = pde_fun(self,x,t,uu,DuDx)
            gamma = 10;
            d = 200;
            
            rho = 13;
            a = 150;
            b = 100;
            alpha = 1.5;
            K = 0.5;

            u = uu(1);
            v = uu(2);
            c = [1; 1;];
            f = [self.D; self.D*d].*DuDx;
            huv = rho.*u.*v./(1+u+K.*(u.^2));
            
            fuv = a-u-huv;
            guv = alpha*(b-v)-huv;
            
            s = [gamma*fuv;
                 gamma*guv];
       end
       
       function u0 = ic_fun(self,x)
            u0 = [0 ; 0];
       end
       
       function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
            % pR is a constant and the initial concentration of phosphorylated protein at the membrane
            % D - diffusion constant
            % pl and ql correspond to the left boundary conditions (x = 0), and 
            % pr and qr correspond to the right boundary condition (x = R).
            pl = [0; 0;];
            ql = [1; 1;];
            pr = [0; 0;];
            qr = [1; 1;];
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

