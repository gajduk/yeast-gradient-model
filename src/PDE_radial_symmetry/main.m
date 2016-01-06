
time_steps = 10000;
end_time = 5000;
spatial_steps = 100;
model = InDirectNegativeFeedbackTwoProteinForms();
output = model.run(time_steps,end_time,spatial_steps);
output.plot_fraction_steady_state()