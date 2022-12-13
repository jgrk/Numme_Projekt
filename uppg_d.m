%c)

clear; close all; clc;
theta = 0:90;
%plot(theta, arrayfun(@f, theta))

%sekantmetoden

V01 = 5; V02 = 15; trunc = 1;

while abs(trunc) > 10^-3

    trunc = f(V01) * (V01 - V02) / (f(V01) - f(V02));
    r = V01 - trunc;
    V02 = V01;
    V01 = r;
    disp(r)
    

end

function trff=f(V0)

%Givna konstanter och funktioner

d=2.37;
m=0.050;
grad=2;
h=1.85;
g=9.82;
Kx=0.001;
Ky=0.01;
be=1.83;

udot=@(u,v)(-(Kx/m)*u*sqrt(u^2+v^2));
vdot=@(u,v)(-g-(Ky/m)*v*sqrt(u^2+v^2));

dt = 6.103515625000000e-06;
u(1)=V0*cos((grad/360)*2*pi);
v(1)=V0*sin((grad/360)*2*pi);
x(1)=0; 
y(1)=h; 

while x(end) < d
    
    x(end+1)=x(end)+u(end)*dt;
    y(end+1)=y(end)+v(end)*dt;
    
    k1=udot(u(end),v(end));
    l1=vdot(u(end),v(end));
    k2=udot(u(end)+k1*dt/2, v(end)+l1*dt/2);
    l2=vdot(u(end)+k1*dt/2, v(end)+l1*dt/2);
    k3=udot(u(end)+k2*dt/2, v(end)+l2*dt/2);
    l3=vdot(u(end)+k2*dt/2, v(end)+l2*dt/2);
    k4=udot(u(end)+k3*dt, v(end)+l3*dt);
    l4=vdot(u(end)+k3*dt, v(end)+l3*dt);

    u(end+1)=u(end)+dt*(k1+2*k2+2*k3+k4)/6;
    v(end+1)=v(end)+dt*(l1+2*l2+2*l3+l4)/6;

end


%Binärsök
if d < x(1)
    mid=1;
end

if d > x(end)
    mid=length(x);
end

hi = length(x);
lo = 1;
while lo <= hi
    mid = fix( ( hi + lo ) * .5 );
    if d < x(mid) 
        hi = mid - 1;
    
    elseif d > x(mid)
        lo = mid + 1;

    else
        break
    end
end

trff = y(mid) - be;

end