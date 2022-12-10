
%Uppgift a)

clear; close all; clc;


%Givna konstanter och funktioner
be=1.83;
d=2.37;
m=0.026;
V0=13;
h=1.85;
grad=5;
g=9.82;
Kx=0.001;
Ky=0.01;
dt=0.1;
tol=10^-5;

udot=@(u,v)(-(Kx/m)*u*sqrt(u^2+v^2));
vdot=@(u,v)(-g-(Ky/m)*v*sqrt(u^2+v^2));

maxiter=50;

for iter = 1:maxiter
    
    clear u; clear v; clear x ; clear y;
    
    u(1)=V0*cos((grad/360)*2*pi);
    v(1)=V0*sin((grad/360)*2*pi);
    
    x(1)=0; 
    y(1)=h; 

    trff_pnkt(1) = 0;
    
    %rk4
    for t=0:dt:.2
    
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
    

    %tidsintervall
    t=0:dt:.2+dt;



    
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
    
    
    %interpolation
    
    t_p = [t(mid-1) t(mid) t(mid+1)]';
    t_mat = [t_p.^0 t_p.^1 t_p.^2];
    x_p = [x(mid-1) x(mid) x(mid+1)]';
    y_p = [y(mid-1) y(mid) y(mid+1)]';
    
    cx=t_mat\x_p; cy=t_mat\y_p;
    x_func = @(t) ( cx(1) + cx(2) * t + cx(3) * t^2 );
    y_func = @(t) ( cy(1) + cy(2) * t + cy(3) * t^2 );
    
    xp_func = @(t) ( cx(2) + 2*cx(3) * t );
    
    %newton-raphson
    
    t0 = .2;    trunc = 1;
    
    while abs( trunc ) > tol
        
        trunc = ( x_func( t0 ) - 2.37 ) / xp_func( t0 );
        t1 = t0 - trunc;
        t0 = t1;
        
    end
    
    trff_pnkt(end+1) = y_func(t1);

    if abs( trff_pnkt(end) - trff_pnkt(end-1) ) < tol
        
        svar = trff_pnkt(end);
        break
    end
    
    dt = dt / 2;



end

hold on
plot(t,x,t,y)
plot(t,arrayfun(x_func,t),t,arrayfun(y_func,t))


disp("träffpunkt från marken: " + svar )


%%

%b)

clear; close all; clc;
theta = 0:90;
%plot(theta, arrayfun(@f, theta))

%sekantmetoden

gr1 = 81; gr2 = 82; trunc = 1;

while abs(trunc) > 10^-3

    trunc = f(gr1) * (gr1 - gr2) / (f(gr1) - f(gr2));
    r = gr1 - trunc;
    gr2 = gr1;
    gr1 = r;
    disp(r)
    

end



function trff=f(grad)

%Givna konstanter och funktioner

d=2.37;
m=0.026;
V0=13;
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

for t=0:dt:1.5
    
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


