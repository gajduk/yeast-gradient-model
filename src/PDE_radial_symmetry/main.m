close all
spatial_steps = 21;
time_steps = 100;

end_time = 30;
model = RealParam_YeastMM_ste11_ste7_fus3_ste5();
output1 = model.run(time_steps,end_time,spatial_steps);
output1.plot_all(1)
