close all
spatial_steps = 50;
models = {InDirectNegativeFeedbackTwoProteinForms(),
          DirectNegativeFeedbackTwoProteinForms(),
          BasicModelTwoProteinForms()};
end_times = [7200,30,5];

%for i=1:3
%    end_time = end_times(i);
%    time_steps = 500;
%    output = models{i}.run(time_steps,end_time,spatial_steps);
%    output.plot_all()
%end
InDirectNegativeFeedbackTwoProteinForms.sensitivity_analysis_kp()