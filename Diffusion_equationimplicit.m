%du/dt = Dd2u/dx2

%using forward euler equation for dudt and central differencing for d2udx2

%Heat Difussion in one dimensional wire using explicit schemes

%Range of space and time
clear

L=1.; %Length of wire
T=1.;  %Final time

Tstep = 2500; %Time descretisation
dt=T/Tstep; %Time difference

n=50; %spatial descretisation
dx=L/n; 
cond = 1./4.;   %conductivity
b=cond*dt/(dx*dx);  %stability parameter (b=<1)

%initial temperature of the wire: a sinus (initialisation)

for i= 1:n+1
    x(i)=(i-1)*dx;
    u(i,1)=sin(pi*x(i));
end

%temperature at the boundary(T=0)
for k=1:Tstep+1
    u(1,k) = 0.;
    u(n+1,k)=0.;
    time(k)=(k-1)*dt;
end

%implementation of implicit method

for k=1:Tstep %time loop
    for i=2:n  %spaceloop
        u(i,k+1)=u(i,k)+0.5*b*(u(i-1,k)+u(i+1,k)-2.*u(i,k));
    end
end

figure(1)
plot(x,u(:,1),'-',x,u(:,100),'-',x,u(:,600),'-');
title('Temperature within the implicit method')
xlabel('X')
ylabel('T')


figure(2)
mesh(x,time,u')
title('Temperature within the implicit method')
xlabel('X')
ylabel('T')