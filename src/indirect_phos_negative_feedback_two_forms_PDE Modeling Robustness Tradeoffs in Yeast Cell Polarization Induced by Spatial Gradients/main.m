defaults;
global pr_list
global t_list
pr_list = [];
t_list = [];
[r,t,sol,u,steady_state,half_fraction] = run_model(time_steps,end_time,spatial_steps,R,D,p0,ka,kp,A,phos0,phos_max,phos_regulation_rate);
plot_fraction(r,t,sol,ka,A,R);