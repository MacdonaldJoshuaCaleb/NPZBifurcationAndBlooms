function fig = ScenarioPlots(aa,tite)
% takes model parameter alpha (maximum harmful affect and plot titile as
% inupt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uu = paramfun(p,t)
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
sgtitle(tite,'Interpreter','latex','Fontsize',20)



% do an initial run to find the correct initial condtitions in ideal
% scenario
p01 = [.5,10,.3,.29,.85,4,.15,.5,10,.4,.3,1,.7,4,6];
ts = 0:100;
sol1 = paramfun(p01,ts);

% run again with the desired scenario and plot
p02 = [.5,10,.3,.29,aa,4,.15,.5,10,.4,.3,sol1(end,1),sol1(end,2),sol1(end,3),sol1(end,4)];
sol1 = paramfun(p02,ts);
hold on
grid on
plot(ts,sol1(:,1),'b','linewidth',3)
plot(ts,sol1(:,2),'r','linewidth',3)
xticks([25,50,75,100])
% xticklabels({'3 months','6 months','9 months','1 year'})
legend({'$p(\tau)$','$z(\tau)$'},'Interpreter','latex','Location','northeastoutside')
ylim([0 8])

hold off
xlabel(' ')
ylabel('Concentration')
set(gca,'FontSize',20,'XTickLabel',[])

nexttile
hold on
grid on
plot(ts,sol1(:,3),'g','linewidth',3)
plot(ts,sol1(:,4),'k','linewidth',3)
hold off
xticks([25,50,75,100])
% xticklabels({'3 months','6 months','9 months','1 year'})
legend({'$n(\tau)$','$d(\tau)$'},'Interpreter','latex','Location','northeastoutside')
ylim([0 16])
xlabel(' ')
ylabel('Concentration')
set(gca,'FontSize',20,'XTickLabel',[])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the persistence numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
grid on
hold on
plot(ts,sol1(:,5),'b-.','linewidth',3)
plot(ts,sol1(:,6),'r-.','linewidth',3)
yline(1,'k--','linewidth',3)
text(40,1.2,'Uncertain persistence threshold','Fontsize',20,'Interpreter','latex')
legend({'$\mathcal{P}_0^p(\tau)$','$\mathcal{P}_0^z(\tau)$'},'Interpreter','latex','Location','northeastoutside')
ylim([0 5])
xticks([25,50,75,100])
ylabel('Persistence number')
xticklabels({'3 months','6 months','9 months','1 year'})
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
plot(ts,1-sol1(:,6)./sol1(:,5),'Color',	[1,0,1],'linewidth',3);
ylim([-1 1])
ylabel('Balance')
xticks([25,50,75,100])
yticks([-1,0,1])
text(42.5,.8,'Favor phytoplankton','Color','blue','Fontsize',20)
text(42.5,-.8,'Favor zooplankton','Color','red','Fontsize',20)
xticklabels({'3 months','6 months','9 months','1 year'})
xtickangle(45)
set(gca,'FontSize',20)
end
