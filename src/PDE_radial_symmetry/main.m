close all
time_steps = 500;
end_time = 30;
spatial_steps = 100;
model = DirectNegativeFeedbackTwoProteinForms();
output = model.run(time_steps,end_time,spatial_steps);
output.plot_all()