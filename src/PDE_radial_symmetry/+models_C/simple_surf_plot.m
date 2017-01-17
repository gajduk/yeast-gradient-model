function simple_surf_plot()
    plotter = models_C.Model1910Plotter();
    time_steps = 2000;
    spatial_steps = 101;
    end_time = 10;
    model1 = models_C.Model1910_C1_and_C2(2);
    [xmesh,tspan,res] = model1.run(time_steps,end_time,spatial_steps);
    plotter.plot_all_surf(xmesh,tspan,res)
    export_fig('simple_surf_plot.png') 
end