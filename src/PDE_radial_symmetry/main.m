close all
spatial_steps = 50;
end_time = 100;
time_steps = 1000;

model = MMInDirectNegativeFeedbackTwoProteinForms();
output = model.run(time_steps,end_time,spatial_steps);
output.plot_all()