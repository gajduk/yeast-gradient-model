close all
spatial_steps = 101;
time_steps = 1000;

end_time = 30;
model = YeastMM_ste11_ste7_fus3_ste5();
output = model.run(time_steps,end_time,spatial_steps);
output.plot_all()

