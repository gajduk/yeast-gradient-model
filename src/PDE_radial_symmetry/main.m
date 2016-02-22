close all
spatial_steps = 100;
model = ThreeProteinsAutocatalytic();
end_time = 2;
time_steps = 1000;
output = model.run(time_steps,end_time,spatial_steps);
figure;
subplot(1,3,1)
output.plot_fraction_internal()
subplot(1,3,2)
output.plot_phosphatase_internal()
                
subplot(1,3,3)
output.plot_erk_internal()