function res= main()
spatial_steps = 101;
time_steps = 1000;
end_time = 10;
model = models_A_and_B.RealParam_YeastMM_ste11_ste7_fus3_ste5();
res = model.run(time_steps,end_time,spatial_steps);
res.plot_all()
end
