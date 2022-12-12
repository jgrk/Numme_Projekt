clc;clear;close all;

%Uppgift a) första utkast

Kx=.001; Ky=.01;




mittpunkt=-0.02;
m=.026; %[m]
V0=13; %[m/s]
g=9.82; %[N/kg]


x2dot=@(xdot,ydot)-(Kx/m)*xdot*sqrt(xdot^2+ydot^2);
y2dot=@(xdot,ydot)-g-(Ky/m)*ydot*sqrt(xdot^2+ydot^2);

trunc_error=1;
tp_y=1;

while abs(trunc_error)>10^-5
    clear y_val; clear x_val; clear xdot_val; clear ydot_val;

    y_val(1)=0; x_val(1)=0;
    ydot_val(1)=V0*sin((5/360)*2*pi); 
    xdot_val(1)=V0*cos((5/360)*2*pi);
    
    h=.01;
    
    for i=1:h:5
        
        y_val(end+1)=y_val(end)+h*ydot_val(end);
        x_val(end+1)=x_val(end)+h*xdot_val(end);
    
        ydot_val(end+1)=ydot_val(end)+h*y2dot(xdot_val(end), ydot_val(end));
        xdot_val(end+1)=xdot_val(end)+h*x2dot(xdot_val(end), ydot_val(end));
    
       
    end
    
    
    t=1:h:5+h;
    hold on
    %plot(t,x_val,t,y_val)
    
    %Erhållt x och y som en funktion av t
    %Behöver interpolera data
    %Behöver hitta x==2.37
    
    t=[t(19) t(20) t(21)]';
    bx=[x_val(19) x_val(20) x_val(21)]';
    by=[y_val(19) y_val(20) y_val(21)]';
    
    A=[t.^0 t.^1 t.^2];
    
    cx=A\bx;
    cy=A\by;
    
    %Index att interpolera 19, 20, 21
    
    t=1:h:5+h;
    
    x_func=@(t)(cx(1)+cx(2)*t+cx(3)*t^2-2.37);
    y_func=@(t)(cy(1)+cy(2)*t+cy(3)*t^2);
    
    xdot_func=@(t)(cx(2)+2*cx(3)*t);
    ydot_func=@(t)(cy(2)+2*cy(3)*t);
    
    %Newton-Raphson
    
    t0=1;
    trunc=1;
    
    while abs(trunc)>10^-8
        
        trunc=x_func(t0)/xdot_func(t0);
        
        t1=t0-trunc;
        t0=t1;
    
    end
    
    plot(t,arrayfun(x_func,t),"g")
    old_tp_y=tp_y;
    
    tp_y=y_func(t1) %träffpunkt relativt y=0 koordinatsystem
    
    delta_y=tp_y-mittpunkt %Träffpunkt relativt bullseye
    
    trunc_error=old_tp_y-tp_y;
    h=h/2;
end


%%

