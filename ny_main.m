clear; close all; clc;
format long e

Kx=0.001;
Ky=0.01;
m=0.026;
V0=13;
g=9.82;
d=2.37;
delta_y=-0.02;

V=@(u,v)sqrt(u^2+v^2);

udot=@(u,v)(-(Kx/m)*u*V(u,v));
vdot=@(u,v)(-g-(Ky/m)*v*V(u,v));

tp_y(1)=1;
h=0.1;
    
for iter=2:20
    
    clear u; clear v; clear x; clear y;
    %BV
    
    u(1)=V0*cos((5/360)*2*pi);
    v(1)=V0*sin((5/360)*2*pi);
    
    x(1)=0;
    y(1)=0;
    
    for t=0:h:0.2
    
        k1=h*u(end);
        l1=h*v(end);
        k2=h*udot(u(end)+k1, v(end)+l1);
        l2=h*vdot(u(end)+k1, v(end)+l1);
        u(end+1)=u(end)+.5*(k1+k2);
        v(end+1)=v(end)+.5*(l1+l2);
    
        x(end+1)=x(end)+h*u(end);
        y(end+1)=y(end)+h*v(end);
    
    end
    
    t=0:h:0.2+h;
    
    plot(t,x)
    
    %för att hita index vars x är närmst 
    
    idx=round(length(t)/2);
    idx_end=length(t);
    idx_begin=1;
    
    while true
        
        g_idx=idx;
        
        if x(idx)<d
            idx_begin=idx;
            idx=round(.5*(idx+idx_end));
        end
      
        if x(idx)>d
            idx_end=idx;
            idx=round(.5*(idx_begin+idx));
        end
            
        if x(idx)==d
            break
        end
        
        if x(idx)==x(g_idx)
            break
        end
    end
    
    %index att interpolera idx-1 idx idx+1
    
    txy=[t(idx-1) t(idx) t(idx+1)]';
    bx=[x(idx-1) x(idx) x(idx+1)]';
    by=[y(idx-1) y(idx) y(idx+1)]';
    
    A=[txy.^0 txy.^1 txy.^2];
    
    cx=A\bx;
    cy=A\by;
    
    x_func=@(t)(cx(1)+cx(2)*t+cx(3)*t^2);
    y_func=@(t)(cy(1)+cy(2)*t+cy(3)*t^2);
    
    xdot_func=@(t)(cx(2)+2*cx(3)*t);
    ydot_func=@(t)(cy(2)+2*cy(3)*t);
    
    hold on
    plot(t,arrayfun(x_func,t),"g")
    plot(t,x,"o")
    plot(t,arrayfun(y_func,t),"b")
    plot(t,y,"x")
    
    t0=1;
    trunc=1;
    
    while abs(trunc) > 10^-8
        
        trunc=x_func(t0) / xdot_func(t0);
        
        t1=t0 - trunc;
        t0=t1;
    
    end
    
    tp_y(iter)=y_func(t1) - delta_y;

    if abs( tp_y(iter) - tp_y(iter - 1) ) < 10^-5
        break
    end
   
    h = h / 2;

end



