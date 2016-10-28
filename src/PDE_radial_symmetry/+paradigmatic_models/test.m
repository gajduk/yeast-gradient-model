function test()
    phos_kin_sensitivity()
end

function simple_surf_plot()
    plotter = paradigmatic_models.Model1910Plotter();
    time_steps = 2000;
    spatial_steps = 101;
    
    end_time = 10;
    model1 = paradigmatic_models.AAModel1910(1);
    [xmesh,tspan,res] = model1.run(time_steps,end_time,spatial_steps);
    plotter.plot_all_surf(xmesh,tspan,res)

end

function phos_kin_sensitivity()
    figure('Position',[100,100,1500,400]);
    subplot(1,3,1)
    test_parameters({'s_phos'},0)
    subplot(1,3,2)
    test_parameters({'s_k_cat'},0)
    subplot(1,3,3)
    test_parameters({'s_phos','s_k_cat'},0)
    export_fig('aa_no_feedback_just_kin_phos.png')
end
function test_parameters(param_list,export)
    plotter = paradigmatic_models.Model1910Plotter();
    
    time_steps = 2000;
    spatial_steps = 101;
    
    end_time = 10;
    
    ks = -3:2:3;
    results = cell(1,length(ks));
    
    parfor i=1:length(ks)
        k = ks(i);
        coef = 5.^k;
        model1 = paradigmatic_models.AAModel1910(1);
        for q=1:length(param_list)
           param = param_list{q};
           model1.(param) = model1.(param)*coef;
        end
        [xmesh,tspan,res] = model1.run(time_steps,end_time,spatial_steps);
        results{i}.res = res;
        results{i}.xmesh = xmesh;
        results{i}.tspan = tspan;
    end
    for i=1:length(ks)
        k = ks(i);
        coef = 5.^k;
        res = results{i}.res;
        xmesh = results{i}.xmesh;
        tspan = results{i}.tspan;
        plotter.plot_gradient(res,5,xmesh,coef);  
        hold on   
    end
    s_title = strjoin(param_list,' ');
    title(s_title)
    legend('show')
    if export == 1
        export_fig(strcat(s_title,'.png')) 
    end
end

