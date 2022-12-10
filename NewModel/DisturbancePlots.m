function [] = SyntheticData_Identifiability
close all
%% function definitions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the model absent disturbance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uu = model(p,t)
    ff = @(t) 1 + (1/2).*sin((2.*pi.*t)./100);

      k=p(1); c=p(2); eps=p(3); gtil = p(4); alpha=p(5); theta=p(6); psi=p(7); g = p(8); xi = p(9); a = p(10); s = p(11);  
      omega = 10; M_theta = 0;
      
         gg = @(t) (1/omega.^2).*t.*exp(-(t)./omega);
          m = max(gg(t));
          theta_d = @(t) (M_theta./m).*gg(t) + theta;
    u0 = [p(12),p(13),p(14),p(15)]; % initial conditions 
     f = @(t,u) [ff(t).*(u(3)./(k+u(3))).*u(1).*(1-u(1)./c) - (u(1).^2./(1+u(1).^2)).*u(2) - eps.*u(1); ... 
               (gtil./ff(t)).*((u(1).^2./(1+u(1).^2))-a.*u(2)-alpha.*(u(1)./(xi+u(1)))).*u(2); ...
                -1.*((u(3)./(k+u(3))).*u(1).*(1-u(1)./c))+s.*(theta_d(t)-u(3))+psi.*u(4); ...
                (1-g).*(u(1).^2./(1+u(1).^2)).*u(2)+eps.*u(1)-psi.*u(4)];
   
   [~,uu] = ode45(f,t,u0);
   R0p = (theta_d(t)./(.3.*(.5+theta_d(t))));
   phat = c.*(1-1./R0p);
   R0z = (phat.*(xi+phat))./(alpha.*(1+phat.^2))';
   R0p = R0p';
   uu = [uu,R0p,R0z];
   
end

function uu = disturbance(p,t)
    ff = @(t) 1 + (1/2).*sin((2.*pi.*t)./100);

      k=p(1); c=p(2); eps=p(3); gtil = p(4); alpha=p(5); theta=p(6); psi=p(7); g = p(8); xi = p(9); a = p(10); s = p(11);  
      omega = p(16); M_theta = p(17);
      
         gg = @(t) (1/omega.^2).*t.*exp(-(t)./omega);
          m = max(gg(t));
          theta_d = @(t) -(M_theta./m).*gg(t) + theta;
    u0 = [p(12),p(13),p(14),p(15)]; % initial conditions 
     f = @(t,u) [ff(t).*(u(3)./(k+u(3))).*u(1).*(1-u(1)./c) - (u(1).^2./(1+u(1).^2)).*u(2) - eps.*u(1); ... 
               (gtil./ff(t)).*((u(1).^2./(1+u(1).^2))-a.*u(2)-alpha.*(u(1)./(xi+u(1)))).*u(2); ...
                -1.*((u(3)./(k+u(3))).*u(1).*(1-u(1)./c))+s.*(theta_d(t)-u(3))+psi.*u(4); ...
                (1-g).*(u(1).^2./(1+u(1).^2)).*u(2)+eps.*u(1)-psi.*u(4)];
   
   [~,uu] = ode45(f,t,u0);
   R0p = (theta_d(t)./(.3.*(.5+theta_d(t))));
   phat = c.*(1-1./R0p);
   R0z = (phat.*(xi+phat))./(alpha.*(1+phat.^2))';
   R0p = R0p';
   uu = [uu,R0p,R0z];
end
%% Monte-carlo simulations 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do a run of the model to get the undisturbed solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw a disturbance start time, peak, and duration  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = randi([2,101],1,1);
omega_r = exprnd(100);
M_theta_r = 4;

ts = 0:100;

ts_dist = ts(r):ts(r) + 100;

p0 = [0.5000   10.0000    0.3000    0.3000    .8    4.0000    0.1500    0.5000   10.0000    0.4000    0.3000    1.3076    0.8305    2.6022    2.9699];

ts = 0:ts_dist(end);
sol1 = model(p0,ts);

p0_dist = [0.5000   10.0000    0.3000    0.3000    .8    4.0000    0.1500    0.5000   10.0000    0.4000    0.3000    sol1(ts(r),1) sol1(ts(r),2)    sol1(ts(r),3)    sol1(ts(r),4),omega_r,M_theta_r];
dist = disturbance(p0_dist,ts_dist);
%% Post proces figure generation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
fig.Position = [86 1 1835 1081];
tlo = tiledlayout(2,2,'TileSpacing', "compact","Padding",'compact');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the time trajectory of the ecosystem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile
sgtitle('Reoligotrophication event - healthy ecosystem','Interpreter','latex','Fontsize',20)





hold on
grid on
plot(ts_dist,sol1(end-100:end,1),'b','linewidth',3)
plot(ts_dist,sol1(end-100:end,2),'r','linewidth',3)
plot(ts_dist,dist(:,1),'b-.','linewidth',3)
plot(ts_dist,dist(:,2),'r-.','linewidth',3)
xticks(ts_dist(1)+[0,25,50,75,100])
% xticklabels({'3 months','6 months','9 months','1 year'})
legend({'$p(\tau)$','$z(\tau)$','$p(\tau)$ (disturbed)','$z(\tau)$ (disturbed)'},'Interpreter','latex','Location','northeastoutside')
ylim([0 8])
xlim([ts_dist(1),ts_dist(end)])
hold off
xlabel(' ')
ylabel('Concentration')
set(gca,'FontSize',20,'XTickLabel',[])

nexttile
hold on
grid on
plot(ts_dist,sol1(end-100:end,3),'g','linewidth',3)
plot(ts_dist,sol1(end-100:end,4),'k','linewidth',3)
plot(ts_dist,dist(:,3),'g-.','linewidth',3)
plot(ts_dist,dist(:,4),'k-.','linewidth',3)
hold off
xticks(ts_dist(1)+[0,25,50,75,100])
xlim([ts_dist(1),ts_dist(end)])
% xticklabels({'3 months','6 months','9 months','1 year'})
legend({'$n(\tau)$','$d(\tau)$','$n(\tau)$ (disturbed)','$d(\tau)$ (disturbed)'},'Interpreter','latex','Location','northeastoutside')
ylim([0 20])
xlabel(' ')
ylabel('Concentration')
set(gca,'FontSize',20,'XTickLabel',[])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the persistence numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
grid on
hold on
plot(ts_dist,sol1(end-100:end,5),'b-.','linewidth',3)
plot(ts_dist,sol1(end-100:end,6),'r-.','linewidth',3)
plot(ts_dist,dist(:,5),'b:','linewidth',3)
plot(ts_dist,dist(:,6),'r:','linewidth',3)
yline(1,'k--','linewidth',3)
text(ts_dist(1)+30,1.2,'Uncertain persistence threshold','Fontsize',20,'Interpreter','latex')
legend({'$\mathcal{P}_0^p(\tau)$','$\mathcal{P}_0^z(\tau)$','$\mathcal{P}_0^p(\tau)$ (disturbed)','$\mathcal{P}_0^z(\tau)$ (disturbed)'},'Interpreter','latex','Location','northeastoutside')
ylim([0 5])
xticks(ts_dist(1)+[0,25,50,75,100])
xlim([ts_dist(1),ts_dist(end)])
ylabel('Persistence number')
xticklabels({'Start','3 months','6 months','9 months','1 year'})
xtickangle(45)
set(gca,'FontSize',20)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % the ecosystem balance
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
nexttile
grid on
hold on
plot(ts_dist,1-sol1(end-100:end,6)./sol1(end-100:end,5),'Color',	[1,0,1],'linewidth',3);
plot(ts_dist,1-dist(:,6)./dist(:,5),'Color',	[1,0,1],'Linestyle','-.','linewidth',3);
ylim([-1 1])
ylabel('Balance')
xticks(ts_dist(1)+[0,25,50,75,100])
xlim([ts_dist(1),ts_dist(end)])
yticks([-1,0,1])
text(ts_dist(1)+42.5,.8,'Favor phytoplankton','Color','blue','Fontsize',20)
text(ts_dist(1)+42.5,-.8,'Favor zooplankton','Color','red','Fontsize',20)
xticklabels({'Start','3 months','6 months','9 months','1 year'})
xtickangle(45)
legend({'$\mathcal{B}(\tau)$','$\mathcal{B}(\tau)$ (disturbed)'},'Interpreter','latex','Location','northeastoutside')
set(gca,'FontSize',20)
end
