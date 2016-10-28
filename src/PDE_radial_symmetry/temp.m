function temp()
    ks = -3:1:3;
    time_steps = 2000;
    spatial_steps = 101;
    end_time = 10;
    for i=1:length(ks)
       k = ks(i);
       coef = 5.^k;
       xmesh = linspace(0,1,spatial_steps);
       m = models.MMInDirectNegativeFeedbackTwoProteinForms();
       m.V_phos_max = m.V_phos_max*coef;
       o = m.run(time_steps,end_time,spatial_steps);
       y = fliplr(o.steady_state);
       title_s = sprintf('%3.3f hf=%.3f',coef,utils.get_hf(xmesh,y));
            
       plot(y./y(1),'DisplayName',title_s)
       hold on
       set(gca,'XTick',[0,99])
       xlim([0,99])
       set(gca,'XTickLabel',{'PM','NM'})
    end
    legend('show')
end

