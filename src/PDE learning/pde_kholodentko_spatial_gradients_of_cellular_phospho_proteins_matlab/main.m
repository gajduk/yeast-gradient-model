defaults;
[r,t,sol,u,steady_state,half_fraction] = run_model(time_steps,end_time,spatial_steps,R,D,p0,ka,kp);
plot_fraction(r,t,sol,ka,kp);