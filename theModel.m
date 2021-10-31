function uu = theModel(x,t)
c = x(1);
k = x(2);
ghat = x(3);
a = x(4);
g = x(5);
s = x(6);
n0 = x(7);
u0 = x(8:10);
f = @(t,u) [(u(3)./(c*k+c*u(3)))*u(1)*(c-u(1)) - (u(1).^2*u(2))./(1+u(1).^2);
            ghat*(u(1).^2*u(2))./(1+u(1).^2)-a*ghat*u(2).^1; % change to 1 for linear loss
            -(u(3)./(c*k+c*u(3)))*u(1)*(c-u(1)) + (1-g)*(u(1).^2*u(2))./(1+u(1).^2) + s*(n0-u(3))];
[~,uu] = ode45(f,t,u0); 

%uu = [uu(end,1),uu(end,2),uu(end,3)];
end

