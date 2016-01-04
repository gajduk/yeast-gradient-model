function res = phosphatase_rate(pa,k_cat,E0,K_M,A)
%PHOSPHATASE_RATE - the activity of the phosphotases increases
%proportionally with $u^p(r,t)$ where
%k_cat = 1 [s^-1]  is the catalytic activity of the phosphatases, 
%E_0 = 75 [nM] is the initial concentration of the phosphotases (note that the concentration doesn't change in the model), 
%K_M = 100 [nM] is the  Michaelis constant and
%A is the strength of the feedback (a free parameter).
res = k_cat*((E0*pa)/(K_M+pa))*((1+A*pa*K_M)/(1+pa*K_M));
end

