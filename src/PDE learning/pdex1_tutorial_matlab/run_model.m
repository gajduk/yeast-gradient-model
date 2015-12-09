x = linspace(0,1,1000);
t = linspace(0,2,1000);

m = 0;
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);

u = sol(:,:,1);
    
surf(x,t,u)    
title('Numerical solution computed with 20 mesh points')
xlabel('Distance x')
ylabel('Time t')