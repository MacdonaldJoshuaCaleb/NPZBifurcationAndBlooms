clear

%c = x(1);
%k = x(2);
%ghat = x(3);
%a = x(4);
%g = x(5);
%s = x(6);
%n0 = x(7);
%u0 = x(8:10);
% if ghat < .29 then max grazing < phytoplankton growth, if ghat > .29 
% max grazing > phytoplankton growth 
% value in original article is ghat = .5 --> eta/gamma = 1.72414
% x = [10,.5,.5,.41,.29,.3,12,1,10,1];
% u1 = theModel(x,linspace(0,150,450));
% 
% x(8) = 2.5;
% x(9) = 2.5;
% u2 = theModel(x,linspace(0,150,450));
% t = linspace(0,150,450);
% figure
% box on
% hold on
% plot(t,u1(:,1),'c',t,u1(:,2),'r','linewidth',2)
% plot(t,u1(:,3),'Color','#D8E219','linewidth',2)
% plot(t,u2(:,1),'c-.',t,u2(:,2),'r-.','linewidth',2)
% plot(t,u2(:,3),'Color','#D8E219','LineStyle','-.','linewidth',2)
% xlabel('time')
% ylabel('concentration')
% set(gca,'FontSize',16)
% hold off
% %legend('phytoplankton','zooplankton','nutrients')
% ylim([0 12])
% 
x = [10,.5,.29,.99,.29,.3,12,1,1,1];
%x = [10,.5,.1,.37,.29,.3,12,1,1,1];
u = theModel(x,linspace(0,300,900));
t2 = linspace(0,300,900);
figure
box on
hold on
plot(t2,u(:,1),'c',t2,u(:,2),'r','linewidth',2)
plot(t2,u(:,3),'Color','#D8E219','linewidth',2)
xlabel('time')
ylabel('concentration')
set(gca,'FontSize',16)
hold off
%legend('phytoplankton','zooplankton','nutrients')
ylim([0 12])
% 
% 
% x = [10,.5,.5,.2,.29,.3,12,1,1,1];
% u = theModel(x,linspace(0,150,450));
% figure
% box on
% hold on
% plot(t,u(:,1),'c',t,u(:,2),'r','linewidth',2)
% plot(t,u(:,3),'Color','#D8E219','linewidth',2)
% xlabel('time')
% ylabel('concentration')
% set(gca,'FontSize',16)
% hold off
% %legend('phytoplankton','zooplankton','nutrients')
% ylim([0 12])
% 
% x = [10,.5,.5,.6,.29,.3,12,1,1,1];
% u = theModel(x,linspace(0,150,450));
% figure
% box on
% hold on
% plot(t,u(:,1),'c',t,u(:,2),'r','linewidth',2)
% plot(t,u(:,3),'Color','#D8E219','linewidth',2)
% xlabel('time')
% ylabel('concentration')
% set(gca,'FontSize',16)
% hold off
% %legend('phytoplankton','zooplankton','nutrients')
% ylim([0 12])

% x = [10,.5,.29,.6,.29,.3,12,1,1,1];
% u = theModel(x,linspace(0,150,450));
% figure
% hold on
% plot(t,u(:,1),'c',t,u(:,2),'r',t,u(:,3),'y','linewidth',2)
% hold off
% %legend('phytoplankton','zooplankton','nutrients')
% ylim([0 12])

as = linspace(.001,.999,600);
n0s = linspace(.001,12,600);
cs = linspace(.001,10,600);
ghats = linspace(.001,.999,600);
its = 200;
ps = linspace(.375,4,its);
ns = linspace(.375,4,its);
zs = linspace(.1,4,its);
[P Z] = meshgrid(ps,zs);
PP = P;
ZZ = Z;
ZZZ = Z;
PPP = P;
for j = 1:its
    for k = 1:its
        %x = [10,.5,.1,.37,.29,.3,12,P(j,k),Z(j,k),1]; % limit cycle quad/zoo dom linear
        %x = [10,.5,.5,.41,.29,.3,12,P(j,k),Z(j,k),1]; % bi-stability
        %x = [10,.5,.5,.2,.29,.3,12,P(j,k),Z(j,k),1]; % zoo dominates -
        %x = [10,.5,.5,.7,.29,.3,12,P(j,k),Z(j,k),1]; % hobf biru linear 
        %done
        x = [10,.5,.29,.999,.29,.3,12,P(j,k),Z(j,k),1]; % hobf biru linear 
        %x = [10,.5,.29,.6,.29,.3,12,P(j,k),Z(j,k),1]; % phyto dominates 
        uu = theModel(x,linspace(0,300,600));
        psol(j,k) = uu(end,1);
        maxp(j,k) = max(uu(:,1));
        HAB(j,k) = psol(j,k)-maxp(j,k);
        if maxp(j,k) <= 3*P(j,k)
            maxp(j,k) = NaN;
            PP(j,k) = NaN;
            ZZ(j,k) = NaN;
            HAB(j,k) = NaN;
            ZZZ(j,k) = NaN;
            PPP(j,k) = NaN;
        end
        if psol(j,k)<= 2*P(j,k)
            HAB(j,k) = 0;
%             ZZZ(j,k) = NaN;
%             PPP(j,k) = NaN;

        end

        if HAB(j,k) < 0
            HAB(j,k) = 1;
        end
        

    end
    
end

figure
contourf(PPP,ZZZ,HAB)
title('HAB')

 figure
 hold on
 contourf(P,Z,psol)
 %plot(.1:.1:5,.1:.1:5,'r--','linewidth',2)
 hold off
 xlabel('Intial Phytoplankton')
 ylabel('Initial Zooplankton')
 set(gca,'FontSize',16)
 
 
figure
hold on
contourf(PP,ZZ,maxp)
%contour(P,Z,psol,'LineColor','k')
%contour(PPP,ZZZ,HAB,'LineColor','k')
hold off
colormap(flipud(autumn))
%colormap(summer)
%cb = colorbar;
xlim([.375, 4])
ylim([.1, 4])

%title('Bloom Occurance, region of bi-stability')
%cb.Label.String = 'Peak plankton abundance';
 xlabel('Intial Phytoplankton')
 ylabel('Initial Zooplankton')
 set(gca,'FontSize',16)
    
% % lower = 0;
% % upper = 0;
% % 
% % Pequib = zeros(600,1);
% % Zequib = zeros(600,1);
% % Nequib = zeros(600,1);
% %      for k = 1:600
% %          x = [10,.5,.2,as(k),.29,.3,12,.1,.1,1];
% %          tmp = theModel(x,linspace(0,300,300));
% %          
% %          x(8) = 2.5;
% %          x(9) = 2.5;
% %         tmp2 = theModel(x,linspace(0,300,300));
% %         bound = abs(tmp2-tmp);
% %         if max(bound) > 0.05
% %              %lb = [k-1,as(k-1)]
% %             break
% %         
% %         end
% %       Pequib(k) = tmp(1);
% %          Zequib(k) = tmp(2);
% %          Nequib(k) = tmp(3);
% % 
% %      end
% %      for k = 600:-1:1
% %          x = [10,.5,.2,as(k),.29,.3,12,.1,.1,1];
% %          tmp = theModel(x,linspace(0,300,300));
% %          
% %          x(8) = 2.5;
% %          x(9) = 2.5;
% %         tmp2 = theModel(x,linspace(0,300,300));
% %         bound = abs(tmp2-tmp);
% %         if max(bound) > 0.05
% %               %ub = [k+1,as(k+1)]
% %             break
% %           
% %         end
%          Pequib(k) = tmp(1);
%          Zequib(k) = tmp(2);
%          Nequib(k) = tmp(3);
%          
% 
%      end
%  for j = 1:600
%      if  Nequib(j) == 0
%          Pequib(j) = NaN;
%          Zequib(j) = NaN;
%          Nequib(j) = NaN;
%      end
%  end
% figure
% hold on
% plot(as(1:length(Pequib)),Pequib,'b',as(1:length(Pequib)),Zequib,'r',as(1:length(Pequib)),Nequib,'g')
% % if lb(1) >0 && ub(1) >0
% %     plot(as(lb(1):ub(1)),Pequib
% % end
% % need another loop for bistability 
% %         for m = 1:600
% %         x = [10,.5,.2,.38,.29,.3,12,xs,1,1];