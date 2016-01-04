defaults;
global pr_list
global t_list
pr_list = [];
t_list = [];
[r,t,sol,u,steady_state,half_fraction] = run_model(time_steps,end_time,spatial_steps,R,D,p0,ka,kp,k_cat,E0,K_M,A);
plot_fraction(r,t,sol,ka,A);
figure;
plot(t_list,pr_list);