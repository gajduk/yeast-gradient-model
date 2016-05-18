
classdef MMInDirectNegativeFeedbackTwoProteinForms < ModelCore
    %InDirectNegativeFeedbackTwoProteinForms
    %with michaelis Menten kinetics
    %See https://www.overleaf.com/2061351vhqrxw#/5221863/ - label michaelis
    %menten first TAC
    properties
        R = 10^-5;%sphere radius
        D = 5*10^-12;%diffusion constant of the phosphorylated protein
        p0 = 300;%[nM] total concentration of protein
        phos0 = 200;%[nM] total concentration of phospatases
        V_kin_max = 50;%[nM/s]
        V_phos_max = 15;%[nM/s]
        V_act_max = 2;%[nM/s]
        V_dea_max = 5;%[nM/s]
        A1 = 50;
        A2 = 2;
        KM1
        KM2
        KM3
        KM4
        KM5
        KM6
    end
        
    methods
        
       function self = MMInDirectNegativeFeedbackTwoProteinForms()
            self.KM1 = self.p0/2;
            self.KM2 = self.p0/2;
            self.KM3 = self.phos0/2;
            self.KM4 = self.phos0/2;
            self.KM5 = self.p0/2;
            self.KM6 = self.phos0/2;
       end
       
       function [c,f,s] = pde_fun(self,x,t,u,DuDx)

            p = u(1);%unphosphorylated form of the protein
            pa = u(2);%active - phosphorylated form of the protein
            phos = u(3);%inactive phophatases
            phos_a = u(4);%active phoshatases
            
            v_kin = self.V_kin_max*(p./(self.KM1+p));
            v_phos = self.V_phos_max*(pa./(self.KM2+pa)).*((1+self.A1*phos_a/self.KM3)./(1+phos_a/self.KM3));

            v_kin = v_kin*self.get_kinease_activity_normalizer(x,self.R);
            
            v_act = self.V_act_max*(phos./(self.KM4+phos)).*((1+self.A2*pa/self.KM5)./(1+pa/self.KM5));
            
            v_dea = self.V_dea_max*(phos_a./(self.KM6+phos_a));

            c = [1; 1; 1; 1];
            f = [self.D; self.D; self.D; self.D].*DuDx;
            s = [2*self.D/x; 2*self.D/x; 2*self.D/x; 2*self.D/x].*DuDx;
            s = s+[ -v_kin+v_phos;%non-active protein
                    +v_kin-v_phos;%active protein
                    -v_act+v_dea;
                    +v_act-v_dea];%phosphotases
       end
       
       function u0 = ic_fun(self,x)
            u0 = [self.p0; 0; self.phos0; 0];
       end
       
       function [pl,ql,pr,qr] = bc_fun(self,xl,ul,xr,ur,t)
            % pR is a constant and the initial concentration of phosphorylated protein at the membrane
            % D - diffusion constant
            % pl and ql correspond to the left boundary conditions (x = 0), and 
            % pr and qr correspond to the right boundary condition (x = R).
            pl = [0; 0; 0; 0];
            ql = [1; 1; 1; 1];
            pr = [0; 0; 0; 0];
            qr = [1; 1; 1; 1];
       end

    end
    methods(Static)
        function sensitivity_analysis_A1()
            close all
            time_steps = 3000;
            end_time = 100;
            spatial_steps = 500;
            A1s = 10.^(-6:.1:6);
            half_fractions = zeros(length(A1s),1);
            highest_fraction = zeros(length(A1s),1);
            transient_time = zeros(length(A1s),1);
            initial_kp = 1;
            
            parfor i=1:length(A1s)
                model = MMInDirectNegativeFeedbackTwoProteinForms();
                model.A1 = initial_kp*A1s(i);
                output = model.run(time_steps,end_time,spatial_steps);
                half_fractions(i) = output.half_fraction;
                highest_fraction(i) = output.highest_fraction();
                transient_time(i) = output.transient_time();
            end
            
            figure('Position', [100, 100, 800, 270]);
            subplot(1,3,1);
            semilogx(A1s,half_fractions);
            xlabel({'Regulation rate strength','relative to default'})
            ylabel({'Half fraction at SS','higher values \Rightarrow steeper gradient'})
            xlim([min(A1s),max(A1s)])
            ax = gca;
            ax.XTick = [10^-6,10^-3,1,10^3,10^6];
            ax.XTickLabel = {'10^{-6}','10^{-3}','10^{0}','10^{3}','10^{6}'};
            grid on
            
            subplot(1,3,2);
            loglog(A1s,highest_fraction);
            xlabel({'Regulation rate strength','relative to default'})
            ylabel({'Peak fraction','Normalized to avg at SS'})
            xlim([min(A1s),max(A1s)])
            ax = gca;
            ax.XTick = [10^-6,10^-3,1,10^3,10^6];
            ax.XTickLabel = {'10^{-6}','10^{-3}','10^{0}','10^{3}','10^{6}'};
            grid on
            
            subplot(1,3,3);
            loglog(A1s,transient_time);
            xlabel({'Regulation rate strength','relative to default'})
            ylabel('Peak time [min]')
            xlim([min(A1s),max(A1s)])
            ax = gca;
            ax.XTick = [10^-6,10^-3,1,10^3,10^6];
            ax.XTickLabel = {'10^{-6}','10^{-3}','10^{0}','10^{3}','10^{6}'};
            grid on
            
            export_fig('sensitivity_analysis_A1','-png','-transparent','-r300')
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

