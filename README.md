# Trapezoidal-Matlab
%%%%%%%%%%%%%%Trapezoid%%%%%%%%%%%%%%%%%%
clear all
h = .002;
T = 0:h:.5;
N = length(T);
Z = zeros(N,1);  
z0 = 2^(1/2);
Z(1) = z0;  
t0 = 0;
tol = 1e-7;
C = zeros(N,1); %count number of Newton iterations done at each time step
for i = 2:N
    
    g0 = z0; %initial guess for Newton iteration. Set to solution computed previously.
    F = g0 - z0 - (h/2)*(my_rhs(t0,g0) + my_rhs(t0-h,z0));
   
    count = 0;
    
    while norm(F) > tol
        Fp = 1+(h/2)*df_dz(t0,g0);
        g0 = g0 - F/Fp;
         F = g0 - z0 - (h/2)*(my_rhs(t0,g0) + my_rhs(t0-h,z0));
       
   
        count = count+1;
    end
    C(i) = count; 
    Z(i) = g0;
    z0 = g0; 
    t0 = t0 +h;
    
end
plot(T,Z,'--p'); hold on
y = (1 + exp((-100*T))).^(1/2);
plot(T,y, 'r')
legend( 'Approx','Exact')
