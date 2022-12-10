function sss = NewNPZ
clear
close all
ps = linspace(.000001,10,50000);
set(0,'defaultTextInterpreter','latex');
% c = 10;
% gamma = .5
% alpha = .5
% a = .1
% theta = 4
% s = .3
% xi = 9
% psi = .1


% zstar = @(p,alpha,xi,a) ((p.^2./(1+p.^2))-((alpha.*p)./(xi+p))).*(1/a);
% dstar = @(p,z,gamma,psi,epsilon) (1./psi).*((1-gamma).*(p.^2./(1+p.^2)).*z + epsilon.*p);
% n_b = @(p,d,c,theta,k,psi,s) (p.*(1-p./c) + s.*(k-theta)-psi.*d);
% nstar = @(p,d,c,theta,k,psi,s) (-n_b(p,d,c,theta,k,psi,s) + sqrt((n_b(p,d,c,theta,k,psi,s)).^2+4.*s.*k.*(s.*theta+psi.*d)))*(1/(2.*s));

zstar = @(n,p,k,c,eps) ((1+p.^2)./p).*((n./(k+n)).*(1-p./c)-eps);
dstar = @(p,z,gamma,psi,epsilon) (1./psi).*((1-gamma).*(p.^2./(1+p.^2)).*z + epsilon.*p);
n_b = @(p,c,theta,k,epsilon,s,gamma) s.*(k-theta) - epsilon.*gamma.*p + gamma.*p.*(1-p./c);
nstar = @(p,c,theta,k,epsilon,s,gamma) (-n_b(p,c,theta,k,epsilon,s,gamma) + sqrt((n_b(p,c,theta,k,epsilon,s,gamma)).^2+4.*s.*k.*(epsilon.*p.*gamma +s.*theta))).*(1./(2.*s));

f = @(p) p.^2./(1+p.^2);
g = @(p,z,a,alpha,xi) a.*z + alpha.*p./(xi+p);
h = @(p,z,a,alpha,xi) f(p) - g(p,z,a,alpha,xi);

function ret = J(n,p,z,param)  
    kk=param(1); c=param(2); epsi=param(3); gtil = param(4); alpha=param(5); theta=param(6); psi=param(7); gg = param(8); xi = param(9); a = param(10); s = param(11);
    ret = [(2*p.^3*z)/(p.^2 + 1)^2 - epsi - (n*(p/c - 1))/(kk + n) - (2*p*z)/(p.^2 + 1) - (n*p)/(c*(kk + n)),-p.^2/(p.^2 + 1),     (n*p*(p/c - 1))/(kk + n).^2 - (p*(p/c - 1))/(kk + n),    0; ...
    gtil*z*((2*p)/(p.^2 + 1) - alpha/(p + xi) - (2*p.^3)/(p.^2 + 1)^2 + (alpha*p)/(p + xi).^2), - gtil*(a*z - p.^2/(p.^2 + 1) + (alpha*p)/(p + xi)) - a*gtil*z, 0,    0; ...
    (n*(p/c - 1))/(kk + n) + (n*p)/(c*(kk + n)), 0, (p*(p/c - 1))/(kk + n) - s - (n*p*(p/c - 1))/(kk + n).^2,  psi; ...
    epsi - (2*p*z*(gg - 1))/(p.^2 + 1) + (2*p.^3*z*(gg - 1))/(p.^2 + 1).^2,-(p.^2*(gg - 1))/(p.^2 + 1), 0 , -psi];
end
% toxicity/hyposixa half saturation xi
eps = linspace(.05,10,5000);


for j = 1:length(eps)
% zs = zstar(ps,1.25,eps(j),.3);
% ds = dstar(ps,zs,.5,.15,.3);
% ns = nstar(ps,ds,10,12,.5,.15,.3);
% hs = h(ns,ps,zs,.5,10,.3);

ns = nstar(ps,10,12,.5,.3,.3,.5);
zs = zstar(ns,ps,.5,10,.3);
hs = h(ps,zs,.4,1.25,eps(j));

% J = @(n,p,z,xi) [(2*p*z)/(p^2 + 1) - (2*p^3*z)/(p^2 + 1)^2 - (n*(p/c - 1))/(k + n) - eps - (n*p)/(c*(k + n)),p^2/(p^2 + 1),(n*p*(p/c - 1))/(k + n)^2 - (p*(p/c - 1))/(k + n),0;
%                  gtil*z*((2*p)/(p^2 + 1) - alpha/(p + xi) - (2*p^3)/(p^2 + 1)^2 + (alpha*p)/(p + xi)^2), -gtil*(a*z - p^2/(p^2 + 1) + (alpha*p)/(p + xi)) - a*gtil*z,0,0;
%                  (n*(p/c - 1))/(k + n) + (n*p)/(c*(k + n)),0,(p*(p/c - 1))/(k + n) - s - (n*p*(p/c - 1))/(k + n)^2,psi;
%                  eps - (2*p*z*(g - 1))/(p^2 + 1) + (2*p^3*z*(g - 1))/(p^2 + 1)^2,-(p^2*(g - 1))/(p^2 + 1),0,-psi]
                 
  
    


count = 1;
tempp = [0,0,0];
tempz = [0,0,0];
tempd = [0,0,0];
tempn = [0,0,0];
jdelim = 0;
for k = 1:length(ps)-1
    if hs(k) > 0 && hs(k+1) < 0 && count == 1
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),1.0,eps(j),.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
%         tempn(count) = nstar(ps(k),tempd(count),10,12,.5,.15,.3);
        tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
        count = count+1;
    elseif hs(k) < 0 && hs(k+1) > 0 && count == 1
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),1.0,eps(j),.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
%         tempn(count) = nstar(ps(k),tempd(count),10,12,.5,.15,.3);
        tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
        count = count+1;
    elseif hs(k) > 0 && hs(k+1) < 0 && count == 2
            tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
        count = count+1;
    elseif hs(k) < 0 && hs(k+1) > 0 && count == 2        
        tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
        count = count+1;
    elseif hs(k) > 0 && hs(k+1) < 0 && count == 3     
        tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
    elseif hs(k) < 0 && hs(k+1) > 0 && count == 3        
        tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
    end
end
%[tempp Ip] = sort(tempp);
tempp(tempp == 0) = NaN;
pp(:,j) = tempp;

tempz(tempz == 0) = NaN;
zz(:,j) = tempz;

tempd(tempd == 0) = NaN;
dd(:,j) = tempd;
tempn(tempn == 0) = NaN;
nn(:,j) = tempn;

end
ii = find(pp(1,:) > 0);
if length(ii) > 0
    for j = ii+1:length(eps)
        temp1= pp(:,j);
        [temp1 Ipp] = sort(temp1);
        pp(:,j) = temp1;
        temp2 = zz(:,j);
        temp2 = temp2(Ipp);
        zz(:,j) = temp2;
        temp3 = dd(:,j);
        temp3 = temp3(Ipp);
        dd(:,j) = temp3;
        temp4 = nn(:,j);
        temp4 = temp4(Ipp);
        nn(:,j) = temp4;
        
    end
end

idxz1 = find(zz(1,:) <0);
idxz3 = find(zz(3,:) < 0);

idp = find(pp(2,:) > 0);
if length(idp) > 0
pp(:,1:idp(1)-1) = flipud(pp(:,1:idp(1)-1));
nn(:,1:idp(1)-1) = flipud(nn(:,1:idp(1)-1));
zz(:,1:idp(1)-1) = flipud(zz(:,1:idp(1)-1));
dd(:,1:idp(1)-1) = flipud(dd(:,1:idp(1)-1));
end
idx1 = find(pp(1,:) > 0);
idx2 = find(pp(2,:) > 0);
idx3 = find(pp(3,:) > 0);

nnStab = zeros(size(nn));
ppStab = zeros(size(nn));
zzStab = zeros(size(nn));
ddStab = zeros(size(nn));

nnUStab = zeros(size(nn));
ppUStab = zeros(size(nn));
zzUStab = zeros(size(nn));
ddUStab = zeros(size(nn));
for zzz = 1:length(idx1)
    %kk=param(1); c=param(2); epsi=param(3); gtil = param(4); alpha=param(5); theta=param(6); psi=param(7); gg = param(8); xi = param(9); a = param(10); s = param(11);
    param = [.5,10,.3,.29,1.25,12,.15,.5,eps(idx1(zzz)),.4,.3];
    %[.5,10,.3,.29,1.25,12,.15,.5,5.0,.4,.3,1,1,4,0]
    jac = J(nn(1,idx1(zzz)),pp(1,idx1(zzz)),zz(1,idx1(zzz)),param);
    stab = real(eig(jac)); 
    if max(stab) < 0
        nnStab(1,idx1(zzz)) = nn(1,idx1(zzz));
        ppStab(1,idx1(zzz)) = pp(1,idx1(zzz));
        zzStab(1,idx1(zzz)) = zz(1,idx1(zzz));
        ddStab(1,idx1(zzz)) = dd(1,idx1(zzz));
    end
    idxx = find(stab > 0);
    if length(idxx) > 0
        nnUStab(1,idx1(zzz)) = nn(1,idx1(zzz));
        ppUStab(1,idx1(zzz)) = pp(1,idx1(zzz));
        zzUStab(1,idx1(zzz)) = zz(1,idx1(zzz));
        ddUStab(1,idx1(zzz)) = dd(1,idx1(zzz));
    end
end
for zzz = 1:length(idx2)
    param = [.5,10,.3,.29,1.25,12,.15,.5,eps(idx2(zzz)),.4,.3];
    jac = J(nn(2,idx2(zzz)),pp(2,idx2(zzz)),zz(2,idx2(zzz)),param);
    stab = max(real(eig(jac))); 
    if max(stab) < 0
        nnStab(2,idx2(zzz)) = nn(2,idx2(zzz));
        ppStab(2,idx2(zzz)) = pp(2,idx2(zzz));
        zzStab(2,idx2(zzz)) = zz(2,idx2(zzz));
        ddStab(2,idx2(zzz)) = dd(2,idx2(zzz));
    end
     idxx = find(stab > 0);
    if length(idxx) > 0
        nnUStab(2,idx2(zzz)) = nn(2,idx2(zzz));
        ppUStab(2,idx2(zzz)) = pp(2,idx2(zzz));
        zzUStab(2,idx2(zzz)) = zz(2,idx2(zzz));
        ddUStab(2,idx2(zzz)) = dd(2,idx2(zzz));
    end
    
end
for zzz = 1:length(idx3)
    param = [.5,10,.3,.29,1.25,12,.15,.5,eps(idx3(zzz)),.4,.3];
    jac = J(nn(3,idx3(zzz)),pp(3,idx3(zzz)),zz(3,idx3(zzz)),param);
    stab = max(real(eig(jac))); 
    if max(stab) < 0
        nnStab(3,idx3(zzz)) = nn(3,idx3(zzz));
        ppStab(3,idx3(zzz)) = pp(3,idx3(zzz));
        zzStab(3,idx3(zzz)) = zz(3,idx3(zzz));
        ddStab(3,idx3(zzz)) = dd(3,idx3(zzz));
    end
    idxx = find(stab > 0);
    if length(idxx) > 0
        nnUStab(3,idx3(zzz)) = nn(3,idx3(zzz));
        ppUStab(3,idx3(zzz)) = pp(3,idx3(zzz));
        zzUStab(3,idx3(zzz)) = zz(3,idx3(zzz));
        ddUStab(3,idx3(zzz)) = dd(3,idx3(zzz));
    end
   
end

R0p = 12./(.3.*(.5+12));
for j = 1:length(eps)
    param = [.5,10,.3,.29,1.25,12,.15,.5,eps(j),.4,.3];
    bcheck(j) = max(real(eig(J(12,10.*(1-.3.*(12+.5)./12),0,param))));
   
    R0z(j) = (10.*(1-1/R0p)*(eps(j)+10.*(1-1/R0p)))./(1.25.*(1+10.^2.*(1-1/R0p)^2));
end
sss = [bcheck;R0z];
color = '#d9e314';
figure
hold on
id1= find(ppStab(1,:) > 0);
id2 = find(ppStab(2,:) > 0);
id3 = find(ppStab(3,:) > 0);
id1U= find(ppUStab(1,:) > 0 & zzUStab(1,:)>0);
id2U = find(ppUStab(2,:) > 0 & zzUStab(2,:) > 0);
id3U = find(ppUStab(3,:) > 0 & zzUStab(3,:) > 0);
bUSidx = find(bcheck > 0);
bSidx = find(bcheck < 0);
plot(eps(id1),ppStab(1,id1),'c','linewidth',2)
plot(eps(id2),ppStab(2,id2),'c','linewidth',2)
plot(eps(id3),ppStab(3,id3),'c','linewidth',2)
plot(eps(id1),zzStab(1,id1),'r','linewidth',2)
plot(eps(id2),zzStab(2,id2),'r','linewidth',2)
plot(eps(id3),zzStab(3,id3),'r','linewidth',2)
plot(eps(id1),ddStab(1,id1),'k','linewidth',2)
plot(eps(id2),ddStab(2,id2),'k','linewidth',2)
plot(eps(id3),ddStab(3,id3),'k','linewidth',2)
plot(eps(id1),nnStab(1,id1),'Color',color,'linewidth',2)
plot(eps(id2),nnStab(2,id2),'Color',color,'linewidth',2)
plot(eps(id3),nnStab(3,id3),'Color',color,'linewidth',2)


plot(eps(id1U),ppUStab(1,id1U),'c--','linewidth',2)
plot(eps(id2U),ppUStab(2,id2U),'c--','linewidth',2)
plot(eps(id3U),ppUStab(3,id3U),'c--','linewidth',2)

plot(eps(id1U),zzUStab(1,id1U),'r--','linewidth',2)
plot(eps(id2U),zzUStab(2,id2U),'r--','linewidth',2)
plot(eps(id3U),zzUStab(3,id3U),'r--','linewidth',2)

plot(eps(id1U),ddUStab(1,id1U),'k--','linewidth',2)
plot(eps(id2U),ddUStab(2,id2U),'k--','linewidth',2)
plot(eps(id3U),ddUStab(3,id3U),'k--','linewidth',2)

plot(eps(id1U),nnUStab(1,id1U),'Color',color,'Linestyle','--','linewidth',2)
plot(eps(id2U),nnUStab(2,id2U),'Color',color,'Linestyle','--','linewidth',2)
plot(eps(id3U),nnUStab(3,id3U),'Color',color,'Linestyle','--','linewidth',2)
plot(eps(bUSidx),ones(1,length(eps(bUSidx))).*10.*(1-.3.*(12+.5)./12),'c--','linewidth',2)
plot(eps(bUSidx),zeros(1,length(eps(bUSidx))),'r--','linewidth',2)
plot(eps(bUSidx),ones(1,length(eps(bUSidx))).*12,'Color',color,'Linestyle','--','linewidth',2)
plot(eps(bUSidx),ones(1,length(eps(bUSidx))).*10.*(1-.3.*(12+.5)./12).*2,'k--','linewidth',2)

plot(eps(bSidx),ones(1,length(eps(bSidx))).*10.*(1-.3.*(12+.5)./12),'c','linewidth',2)
plot(eps(bSidx),zeros(1,length(eps(bSidx))),'r','linewidth',2)
plot(eps(bSidx),ones(1,length(eps(bSidx))).*12,'Color',color,'linewidth',2)
plot(eps(bSidx),ones(1,length(eps(bSidx))).*10.*(1-.3.*(12+.5)./12).*2,'k','linewidth',2)
% plot(eps(id1),zzStab(1,id1),'r','linewidth',2)
% plot(eps(id2),zzStab(2,id2),'r','linewidth',2)
% plot(eps(id3),zzStab(3,id3),'r','linewidth',2)
% plot(eps(id1),ddStab(1,id1),'k','linewidth',2)
% plot(eps(id2),ddStab(2,id2),'k','linewidth',2)
% plot(eps(id3),ddStab(3,id3),'k','linewidth',2)
% plot(eps(id1),nnStab(1,id1),'Color',color,'linewidth',2)
% plot(eps(id2),nnStab(2,id2),'Color',color,'linewidth',2)
% plot(eps(id3),nnStabsignficantly(3,id3),'Color',color,'linewidth',2)
ii2 = find(pp(2,:)>0);
if length(ii2) > 0
    xx = [eps(ii2(1)-1),eps(ii2(end)+1)];
    x = [xx,fliplr(xx)];
    ylb = [0,0];
    yub = [14,14];
    inBetween = [ylb,fliplr(yub)];
    hh = fill(x, inBetween,'k','LineStyle','-','Marker','none');
    set(hh,'facealpha',.2)
    inBetween = [ylb,fliplr(yub)];
end
    xx = [min(eps(bSidx(1)),eps(bSidx(end))),max(eps(bSidx(1)),eps(bSidx(end)))];
    x = [xx,fliplr(xx)];
    ylb = [0, 0];
    yub = [14,14];
    inBetween = [ylb,fliplr(yub)];
    hh = fill(x, inBetween,'k','LineStyle','-','Marker','none');
    set(hh,'facealpha',.2)
    inBetween = [ylb,fliplr(yub)];
% if length(idxz1) > 0
% xlim([eps(idxz1(end)+1) eps(end)])
% end
% if length(idxz3)>0
% xlim([eps(idxz3(end)+1) eps(end)])
% end
ylim([0 inf])
xlim([eps(1) eps(end)])
xlabel('$\xi$, half-sat. harmful affect','fontsize',18)
hold off
set(gca,'FontSize',18)
ss = ppUStab;


figure
hold on
id1= find(ppStab(1,:) > 0);
id2 = find(ppStab(2,:) > 0);
id3 = find(ppStab(3,:) > 0);
id1U= find(ppUStab(1,:) > 0 & zzUStab(1,:)>0);
id2U = find(ppUStab(2,:) > 0 & zzUStab(2,:) > 0);
id3U = find(ppUStab(3,:) > 0 & zzUStab(3,:) > 0);
bUSidx = find(bcheck > 0);
bSidx = find(bcheck < 0);
plot(R0z(id1),ppStab(1,id1),'c','linewidth',2)
plot(R0z(id2),ppStab(2,id2),'c','linewidth',2)
plot(R0z(id3),ppStab(3,id3),'c','linewidth',2)
plot(R0z(id1),zzStab(1,id1),'r','linewidth',2)
plot(R0z(id2),zzStab(2,id2),'r','linewidth',2)
plot(R0z(id3),zzStab(3,id3),'r','linewidth',2)
plot(R0z(id1),ddStab(1,id1),'k','linewidth',2)
plot(R0z(id2),ddStab(2,id2),'k','linewidth',2)
plot(R0z(id3),ddStab(3,id3),'k','linewidth',2)
plot(R0z(id1),nnStab(1,id1),'Color',color,'linewidth',2)
plot(R0z(id2),nnStab(2,id2),'Color',color,'linewidth',2)
plot(R0z(id3),nnStab(3,id3),'Color',color,'linewidth',2)


plot(R0z(id1U),ppUStab(1,id1U),'c--','linewidth',2)
plot(R0z(id2U),ppUStab(2,id2U),'c--','linewidth',2)
plot(R0z(id3U),ppUStab(3,id3U),'c--','linewidth',2)

plot(R0z(id1U),zzUStab(1,id1U),'r--','linewidth',2)
plot(R0z(id2U),zzUStab(2,id2U),'r--','linewidth',2)
plot(R0z(id3U),zzUStab(3,id3U),'r--','linewidth',2)

plot(R0z(id1U),ddUStab(1,id1U),'k--','linewidth',2)
plot(R0z(id2U),ddUStab(2,id2U),'k--','linewidth',2)
plot(R0z(id3U),ddUStab(3,id3U),'k--','linewidth',2)

plot(R0z(id1U),nnUStab(1,id1U),'Color',color,'Linestyle','--','linewidth',2)
plot(R0z(id2U),nnUStab(2,id2U),'Color',color,'Linestyle','--','linewidth',2)
plot(R0z(id3U),nnUStab(3,id3U),'Color',color,'Linestyle','--','linewidth',2)
plot(R0z(bUSidx),ones(1,length(R0z(bUSidx))).*10.*(1-.3.*(12+.5)./12),'c--','linewidth',2)
plot(R0z(bUSidx),zeros(1,length(R0z(bUSidx))),'r--','linewidth',2)
plot(R0z(bUSidx),ones(1,length(R0z(bUSidx))).*12,'Color',color,'Linestyle','--','linewidth',2)
plot(R0z(bUSidx),ones(1,length(R0z(bUSidx))).*10.*(1-.3.*(12+.5)./12).*2,'k--','linewidth',2)

plot(R0z(bSidx),ones(1,length(R0z(bSidx))).*10.*(1-.3.*(12+.5)./12),'c','linewidth',2)
plot(R0z(bSidx),zeros(1,length(R0z(bSidx))),'r','linewidth',2)
plot(R0z(bSidx),ones(1,length(R0z(bSidx))).*12,'Color',color,'linewidth',2)
plot(R0z(bSidx),ones(1,length(R0z(bSidx))).*10.*(1-.3.*(12+.5)./12).*2,'k','linewidth',2)
% plot(eps(id1),zzStab(1,id1),'r','linewidth',2)
% plot(eps(id2),zzStab(2,id2),'r','linewidth',2)
% plot(eps(id3),zzStab(3,id3),'r','linewidth',2)
% plot(eps(id1),ddStab(1,id1),'k','linewidth',2)
% plot(eps(id2),ddStab(2,id2),'k','linewidth',2)
% plot(eps(id3),ddStab(3,id3),'k','linewidth',2)
% plot(eps(id1),nnStab(1,id1),'Color',color,'linewidth',2)
% plot(eps(id2),nnStab(2,id2),'Color',color,'linewidth',2)
% plot(eps(id3),nnStabsignficantly(3,id3),'Color',color,'linewidth',2)
ii2 = find(pp(2,:)>0);
if length(ii2) > 0
    xx = [min(R0z(ii2(1)-1),R0z(ii2(end)+1)),max(R0z(ii2(1)-1),R0z(ii2(end)+1))];
    x = [xx,fliplr(xx)];
    ylb = [0,0];
    yub = [14,14];
    inBetween = [ylb,fliplr(yub)];
    hh = fill(x, inBetween,'k','LineStyle','-','Marker','none');
    set(hh,'facealpha',.2)
    inBetween = [ylb,fliplr(yub)];
end
   xx = [min(R0z(bSidx(1)),R0z(bSidx(end))),max(R0z(bSidx(1)),R0z(bSidx(end)))];
    x = [xx,fliplr(xx)];
    ylb = [0, 0];
    yub = [14,14];
    inBetween = [ylb,fliplr(yub)];
    hh = fill(x, inBetween,'k','LineStyle','-','Marker','none');
    set(hh,'facealpha',.2)
    inBetween = [ylb,fliplr(yub)];
% if length(idxz1) > 0
% xlim([eps(idxz1(end)+1) eps(end)])
% end
% if length(idxz3)>0
% xlim([eps(idxz3(end)+1) eps(end)])
% end
ylim([0 inf])
xlim([min(R0z),max(R0z)])
xlabel('$\mathcal{R}_0^{z}$ as function of $\xi$','fontsize',18)
set(gca,'FontSize',18)
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps = linspace(0.05,3,5000);


for j = 1:length(eps)
% zs = zstar(ps,1.25,eps(j),.3);
% ds = dstar(ps,zs,.5,.15,.3);
% ns = nstar(ps,ds,10,12,.5,.15,.3);
% hs = h(ns,ps,zs,.5,10,.3);

ns = nstar(ps,10,12,.5,.3,.3,.5);
zs = zstar(ns,ps,.5,10,.3);
hs = h(ps,zs,.4,eps(j),7);

% J = @(n,p,z,xi) [(2*p*z)/(p^2 + 1) - (2*p^3*z)/(p^2 + 1)^2 - (n*(p/c - 1))/(k + n) - eps - (n*p)/(c*(k + n)),p^2/(p^2 + 1),(n*p*(p/c - 1))/(k + n)^2 - (p*(p/c - 1))/(k + n),0;
%                  gtil*z*((2*p)/(p^2 + 1) - alpha/(p + xi) - (2*p^3)/(p^2 + 1)^2 + (alpha*p)/(p + xi)^2), -gtil*(a*z - p^2/(p^2 + 1) + (alpha*p)/(p + xi)) - a*gtil*z,0,0;
%                  (n*(p/c - 1))/(k + n) + (n*p)/(c*(k + n)),0,(p*(p/c - 1))/(k + n) - s - (n*p*(p/c - 1))/(k + n)^2,psi;
%                  eps - (2*p*z*(g - 1))/(p^2 + 1) + (2*p^3*z*(g - 1))/(p^2 + 1)^2,-(p^2*(g - 1))/(p^2 + 1),0,-psi]
                 
  
    


count = 1;
tempp = [0,0,0];
tempz = [0,0,0];
tempd = [0,0,0];
tempn = [0,0,0];
jdelim = 0;
for k = 1:length(ps)-1
    if hs(k) > 0 && hs(k+1) < 0 && count == 1
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),1.0,eps(j),.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
%         tempn(count) = nstar(ps(k),tempd(count),10,12,.5,.15,.3);
        tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
        count = count+1;
    elseif hs(k) < 0 && hs(k+1) > 0 && count == 1
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),1.0,eps(j),.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
%         tempn(count) = nstar(ps(k),tempd(count),10,12,.5,.15,.3);
        tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
        count = count+1;
    elseif hs(k) > 0 && hs(k+1) < 0 && count == 2
            tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
        count = count+1;
    elseif hs(k) < 0 && hs(k+1) > 0 && count == 2        
        tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
        count = count+1;
    elseif hs(k) > 0 && hs(k+1) < 0 && count == 3     
        tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
    elseif hs(k) < 0 && hs(k+1) > 0 && count == 3        
        tempp(count) = ps(k);
        tempn(count) = nstar(ps(k),10,12,.5,.3,.3,.5);
        tempz(count) = zstar(tempn(count),ps(k),.5,10,.3);
        tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
    end
end
%[tempp Ip] = sort(tempp);
tempp(tempp == 0) = NaN;
pp(:,j) = tempp;

tempz(tempz == 0) = NaN;
zz(:,j) = tempz;

tempd(tempd == 0) = NaN;
dd(:,j) = tempd;
tempn(tempn == 0) = NaN;
nn(:,j) = tempn;

end
ii = find(pp(1,:) > 0);
if length(ii) > 0
    for j = ii+1:length(eps)
        temp1= pp(:,j);
        [temp1 Ipp] = sort(temp1);
        pp(:,j) = temp1;
        temp2 = zz(:,j);
        temp2 = temp2(Ipp);
        zz(:,j) = temp2;
        temp3 = dd(:,j);
        temp3 = temp3(Ipp);
        dd(:,j) = temp3;
        temp4 = nn(:,j);
        temp4 = temp4(Ipp);
        nn(:,j) = temp4;
        
    end
end

idxz1 = find(zz(1,:) <0);
idxz3 = find(zz(3,:) < 0);

idp = find(pp(2,:) > 0);
if length(idp) > 0
pp(:,idp(end)+1:end) = flipud(pp(:,idp(end)+1:end));
nn(:,idp(end)+1:end) = flipud(nn(:,idp(end)+1:end));
zz(:,idp(end)+1:end) = flipud(zz(:,idp(end)+1:end));
dd(:,idp(end)+1:end) = flipud(dd(:,idp(end)+1:end));
end
idx1 = find(pp(1,:) > 0);
idx2 = find(pp(2,:) > 0);
idx3 = find(pp(3,:) > 0);

nnStab = zeros(size(nn));
ppStab = zeros(size(nn));
zzStab = zeros(size(nn));
ddStab = zeros(size(nn));

nnUStab = zeros(size(nn));
ppUStab = zeros(size(nn));
zzUStab = zeros(size(nn));
ddUStab = zeros(size(nn));
for zzz = 1:length(idx1)
    %kk=param(1); c=param(2); epsi=param(3); gtil = param(4); alpha=param(5); theta=param(6); psi=param(7); gg = param(8); xi = param(9); a = param(10); s = param(11);
    param = [.5,10,.3,.29,eps(idx1(zzz)),12,.15,.5,7,.4,.3];
    %[.5,10,.3,.29,1.25,12,.15,.5,5.0,.4,.3,1,1,4,0]
    jac = J(nn(1,idx1(zzz)),pp(1,idx1(zzz)),zz(1,idx1(zzz)),param);
    stab = real(eig(jac)); 
    if max(stab) < 0
        nnStab(1,idx1(zzz)) = nn(1,idx1(zzz));
        ppStab(1,idx1(zzz)) = pp(1,idx1(zzz));
        zzStab(1,idx1(zzz)) = zz(1,idx1(zzz));
        ddStab(1,idx1(zzz)) = dd(1,idx1(zzz));
    end
    idxx = find(stab > 0);
    if length(idxx) > 0
        nnUStab(1,idx1(zzz)) = nn(1,idx1(zzz));
        ppUStab(1,idx1(zzz)) = pp(1,idx1(zzz));
        zzUStab(1,idx1(zzz)) = zz(1,idx1(zzz));
        ddUStab(1,idx1(zzz)) = dd(1,idx1(zzz));
    end
end
for zzz = 1:length(idx2)
    param = [.5,10,.3,.29,eps(idx2(zzz)),12,.15,.5,7,.4,.3];
    jac = J(nn(2,idx2(zzz)),pp(2,idx2(zzz)),zz(2,idx2(zzz)),param);
    stab = max(real(eig(jac))); 
    if max(stab) < 0
        nnStab(2,idx2(zzz)) = nn(2,idx2(zzz));
        ppStab(2,idx2(zzz)) = pp(2,idx2(zzz));
        zzStab(2,idx2(zzz)) = zz(2,idx2(zzz));
        ddStab(2,idx2(zzz)) = dd(2,idx2(zzz));
    end
     idxx = find(stab > 0);
    if length(idxx) > 0
        nnUStab(2,idx2(zzz)) = nn(2,idx2(zzz));
        ppUStab(2,idx2(zzz)) = pp(2,idx2(zzz));
        zzUStab(2,idx2(zzz)) = zz(2,idx2(zzz));
        ddUStab(2,idx2(zzz)) = dd(2,idx2(zzz));
    end
    
end
for zzz = 1:length(idx3)
    param = [.5,10,.3,.29,eps(idx3(zzz)),12,.15,.5,7,.4,.3];
    jac = J(nn(3,idx3(zzz)),pp(3,idx3(zzz)),zz(3,idx3(zzz)),param);
    stab = max(real(eig(jac))); 
    if max(stab) < 0
        nnStab(3,idx3(zzz)) = nn(3,idx3(zzz));
        ppStab(3,idx3(zzz)) = pp(3,idx3(zzz));
        zzStab(3,idx3(zzz)) = zz(3,idx3(zzz));
        ddStab(3,idx3(zzz)) = dd(3,idx3(zzz));
    end
    idxx = find(stab > 0);
    if length(idxx) > 0
        nnUStab(3,idx3(zzz)) = nn(3,idx3(zzz));
        ppUStab(3,idx3(zzz)) = pp(3,idx3(zzz));
        zzUStab(3,idx3(zzz)) = zz(3,idx3(zzz));
        ddUStab(3,idx3(zzz)) = dd(3,idx3(zzz));
    end
   
end

for j = 1:length(eps)
    param = [.5,10,.3,.29,eps(j),12,.15,.5,7,.4,.3];
    bcheck(j) = max(real(eig(J(12,10.*(1-.3.*(12+.5)./12),0,param))));
    R0z(j) = (10.*(1-1/R0p)*(7+10.*(1-1/R0p)))./(eps(j).*(1+10.^2.*(1-1/R0p)^2));
end
        
    
figure
hold on
id1= find(ppStab(1,:) > 0);
id2 = find(ppStab(2,:) > 0);
id3 = find(ppStab(3,:) > 0);
id1U= find(ppUStab(1,:) > 0 & zzUStab(1,:)>0);
id2U = find(ppUStab(2,:) > 0 & zzUStab(2,:) > 0);
id3U = find(ppUStab(3,:) > 0 & zzUStab(3,:) > 0);
bUSidx = find(bcheck > 0);
bSidx = find(bcheck < 0);
plot(eps(id1),ppStab(1,id1),'c','linewidth',2)
plot(eps(id2),ppStab(2,id2),'c','linewidth',2)
plot(eps(id3),ppStab(3,id3),'c','linewidth',2)
plot(eps(id1),zzStab(1,id1),'r','linewidth',2)
plot(eps(id2),zzStab(2,id2),'r','linewidth',2)
plot(eps(id3),zzStab(3,id3),'r','linewidth',2)
plot(eps(id1),ddStab(1,id1),'k','linewidth',2)
plot(eps(id2),ddStab(2,id2),'k','linewidth',2)
plot(eps(id3),ddStab(3,id3),'k','linewidth',2)
plot(eps(id1),nnStab(1,id1),'Color',color,'linewidth',2)
plot(eps(id2),nnStab(2,id2),'Color',color,'linewidth',2)
plot(eps(id3),nnStab(3,id3),'Color',color,'linewidth',2)


plot(eps(id1U),ppUStab(1,id1U),'c--','linewidth',2)
plot(eps(id2U),ppUStab(2,id2U),'c--','linewidth',2)
plot(eps(id3U),ppUStab(3,id3U),'c--','linewidth',2)

plot(eps(id1U),zzUStab(1,id1U),'r--','linewidth',2)
plot(eps(id2U),zzUStab(2,id2U),'r--','linewidth',2)
plot(eps(id3U),zzUStab(3,id3U),'r--','linewidth',2)

plot(eps(id1U),ddUStab(1,id1U),'k--','linewidth',2)
plot(eps(id2U),ddUStab(2,id2U),'k--','linewidth',2)
plot(eps(id3U),ddUStab(3,id3U),'k--','linewidth',2)

plot(eps(id1U),nnUStab(1,id1U),'Color',color,'Linestyle','--','linewidth',2)
plot(eps(id2U),nnUStab(2,id2U),'Color',color,'Linestyle','--','linewidth',2)
plot(eps(id3U),nnUStab(3,id3U),'Color',color,'Linestyle','--','linewidth',2)
plot(eps(bUSidx),ones(1,length(eps(bUSidx))).*10.*(1-.3.*(12+.5)./12),'c--','linewidth',2)
plot(eps(bUSidx),zeros(1,length(eps(bUSidx))),'r--','linewidth',2)
plot(eps(bUSidx),ones(1,length(eps(bUSidx))).*12,'Color',color,'Linestyle','--','linewidth',2)
plot(eps(bUSidx),ones(1,length(eps(bUSidx))).*10.*(1-.3.*(12+.5)./12).*2,'k--','linewidth',2)

plot(eps(bSidx),ones(1,length(eps(bSidx))).*10.*(1-.3.*(12+.5)./12),'c','linewidth',2)
plot(eps(bSidx),zeros(1,length(eps(bSidx))),'r','linewidth',2)
plot(eps(bSidx),ones(1,length(eps(bSidx))).*12,'Color',color,'linewidth',2)
plot(eps(bSidx),ones(1,length(eps(bSidx))).*10.*(1-.3.*(12+.5)./12).*2,'k','linewidth',2)
% plot(eps(id1),zzStab(1,id1),'r','linewidth',2)
% plot(eps(id2),zzStab(2,id2),'r','linewidth',2)
% plot(eps(id3),zzStab(3,id3),'r','linewidth',2)
% plot(eps(id1),ddStab(1,id1),'k','linewidth',2)
% plot(eps(id2),ddStab(2,id2),'k','linewidth',2)
% plot(eps(id3),ddStab(3,id3),'k','linewidth',2)
% plot(eps(id1),nnStab(1,id1),'Color',color,'linewidth',2)
% plot(eps(id2),nnStab(2,id2),'Color',color,'linewidth',2)
% plot(eps(id3),nnStabsignficantly(3,id3),'Color',color,'linewidth',2)
ii2 = find(pp(2,:)>0);
if length(ii2) > 0
    xx = [eps(ii2(1)-1),eps(ii2(end)+1)];
    x = [xx,fliplr(xx)];
    ylb = [0,0];
    yub = [14,14];
    inBetween = [ylb,fliplr(yub)];
    hh = fill(x, inBetween,'k','LineStyle','-','Marker','none');
    set(hh,'facealpha',.2)
    inBetween = [ylb,fliplr(yub)];
end
   xx = [min(eps(bSidx(1)),eps(bSidx(end))),max(eps(bSidx(1)),eps(bSidx(end)))];
    x = [xx,fliplr(xx)];
    ylb = [0, 0];
    yub = [14,14];
    inBetween = [ylb,fliplr(yub)];
    hh = fill(x, inBetween,'k','LineStyle','-','Marker','none');
    set(hh,'facealpha',.2)
    inBetween = [ylb,fliplr(yub)];
% if length(idxz1) > 0
% xlim([eps(idxz1(end)+1) eps(end)])
% end
% if length(idxz3)>0
% xlim([eps(idxz3(end)+1) eps(end)])
% end
ylim([0 inf])
xlabel('$\tilde{\alpha}$, maximum harmful affect','fontsize',18)
set(gca,'FontSize',18)
hold off


figure
hold on
id1= find(ppStab(1,:) > 0);
id2 = find(ppStab(2,:) > 0);
id3 = find(ppStab(3,:) > 0);
id1U= find(ppUStab(1,:) > 0 & zzUStab(1,:)>0);
id2U = find(ppUStab(2,:) > 0 & zzUStab(2,:) > 0);
id3U = find(ppUStab(3,:) > 0 & zzUStab(3,:) > 0);
bUSidx = find(bcheck > 0);
bSidx = find(bcheck < 0);
plot(R0z(id1),ppStab(1,id1),'c','linewidth',2)
plot(R0z(id2),ppStab(2,id2),'c','linewidth',2)
plot(R0z(id3),ppStab(3,id3),'c','linewidth',2)
plot(R0z(id1),zzStab(1,id1),'r','linewidth',2)
plot(R0z(id2),zzStab(2,id2),'r','linewidth',2)
plot(R0z(id3),zzStab(3,id3),'r','linewidth',2)
plot(R0z(id1),ddStab(1,id1),'k','linewidth',2)
plot(R0z(id2),ddStab(2,id2),'k','linewidth',2)
plot(R0z(id3),ddStab(3,id3),'k','linewidth',2)
plot(R0z(id1),nnStab(1,id1),'Color',color,'linewidth',2)
plot(R0z(id2),nnStab(2,id2),'Color',color,'linewidth',2)
plot(R0z(id3),nnStab(3,id3),'Color',color,'linewidth',2)


plot(R0z(id1U),ppUStab(1,id1U),'c--','linewidth',2)
plot(R0z(id2U),ppUStab(2,id2U),'c--','linewidth',2)
plot(R0z(id3U),ppUStab(3,id3U),'c--','linewidth',2)

plot(R0z(id1U),zzUStab(1,id1U),'r--','linewidth',2)
plot(R0z(id2U),zzUStab(2,id2U),'r--','linewidth',2)
plot(R0z(id3U),zzUStab(3,id3U),'r--','linewidth',2)

plot(R0z(id1U),ddUStab(1,id1U),'k--','linewidth',2)
plot(R0z(id2U),ddUStab(2,id2U),'k--','linewidth',2)
plot(R0z(id3U),ddUStab(3,id3U),'k--','linewidth',2)

plot(R0z(id1U),nnUStab(1,id1U),'Color',color,'Linestyle','--','linewidth',2)
plot(R0z(id2U),nnUStab(2,id2U),'Color',color,'Linestyle','--','linewidth',2)
plot(R0z(id3U),nnUStab(3,id3U),'Color',color,'Linestyle','--','linewidth',2)
plot(R0z(bUSidx),ones(1,length(R0z(bUSidx))).*10.*(1-.3.*(12+.5)./12),'c--','linewidth',2)
plot(R0z(bUSidx),zeros(1,length(R0z(bUSidx))),'r--','linewidth',2)
plot(R0z(bUSidx),ones(1,length(R0z(bUSidx))).*12,'Color',color,'Linestyle','--','linewidth',2)
plot(R0z(bUSidx),ones(1,length(R0z(bUSidx))).*10.*(1-.3.*(12+.5)./12).*2,'k--','linewidth',2)

plot(R0z(bSidx),ones(1,length(R0z(bSidx))).*10.*(1-.3.*(12+.5)./12),'c','linewidth',2)
plot(R0z(bSidx),zeros(1,length(R0z(bSidx))),'r','linewidth',2)
plot(R0z(bSidx),ones(1,length(R0z(bSidx))).*12,'Color',color,'linewidth',2)
plot(R0z(bSidx),ones(1,length(R0z(bSidx))).*10.*(1-.3.*(12+.5)./12).*2,'k','linewidth',2)
% plot(eps(id1),zzStab(1,id1),'r','linewidth',2)
% plot(eps(id2),zzStab(2,id2),'r','linewidth',2)
% plot(eps(id3),zzStab(3,id3),'r','linewidth',2)
% plot(eps(id1),ddStab(1,id1),'k','linewidth',2)
% plot(eps(id2),ddStab(2,id2),'k','linewidth',2)
% plot(eps(id3),ddStab(3,id3),'k','linewidth',2)
% plot(eps(id1),nnStab(1,id1),'Color',color,'linewidth',2)
% plot(eps(id2),nnStab(2,id2),'Color',color,'linewidth',2)
% plot(eps(id3),nnStabsignficantly(3,id3),'Color',color,'linewidth',2)
ii2 = find(pp(2,:)>0);
if length(ii2) > 0
    xx = [min(R0z(ii2(1)-1),R0z(ii2(end)+1)),max(R0z(ii2(1)-1),R0z(ii2(end)+1))];
    x = [xx,fliplr(xx)];
    ylb = [0,0];
    yub = [14,14];
    inBetween = [ylb,fliplr(yub)];
    hh = fill(x, inBetween,'k','LineStyle','-','Marker','none');
    set(hh,'facealpha',.2)
    inBetween = [ylb,fliplr(yub)];
end
   xx = [min(R0z(bSidx(1)),R0z(bSidx(end))),max(R0z(bSidx(1)),R0z(bSidx(end)))];
    x = [xx,fliplr(xx)];
    ylb = [0, 0];
    yub = [14,14];
    inBetween = [ylb,fliplr(yub)];
    hh = fill(x, inBetween,'k','LineStyle','-','Marker','none');
    set(hh,'facealpha',.2)
    inBetween = [ylb,fliplr(yub)];
% if length(idxz1) > 0
% xlim([eps(idxz1(end)+1) eps(end)])
% end
% if length(idxz3)>0
% xlim([eps(idxz3(end)+1) eps(end)])
% end
ylim([0 inf])
xlim([min(R0z),3])
xlabel('$\mathcal{R}_0^{z}$','fontsize',18)
set(gca,'FontSize',18)
hold off



% clear all
% set(0,'defaultTextInterpreter','latex');
% ps = linspace(.000001,10,50000);
% 
% zstar = @(p,alpha,xi,a) ((p.^2./(1+p.^2))-((alpha.*p)./(xi+p))).*(1/a);
% dstar = @(p,z,gamma,psi,epsilon) (1./psi).*((1-gamma).*(p.^2./(1+p.^2)).*z + epsilon.*p);
% n_b = @(p,d,c,theta,k,psi,s) (p.*(1-p./c) + s.*(k-theta)-psi.*d);
% nstar = @(p,d,c,theta,k,psi,s) (-n_b(p,d,c,theta,k,psi,s) + sqrt((n_b(p,d,c,theta,k,psi,s)).^2+4.*s.*k.*(s.*theta+psi.*d)))*(1/(2.*s));
% f = @(n,p,k,c) (n./(k+n)).*p.*(1-p./c);
% g = @(p,z,epsilon) (p.^2./(1+p.^2)).*z + epsilon.*p;
% h = @(n,p,z,k,c,epsilon) f(n,p,k,c) - g(p,z,epsilon);
% 
% 
% % toxicity/hyposixa half saturation 
% eps = linspace(.001,2,1000);
% for j = 1:length(eps)
% zs = zstar(ps,eps(j),4.5,.4);
% ds = dstar(ps,zs,.5,.15,.3);
% ns = nstar(ps,ds,10,4,.5,.15,.3);
% hs = h(ns,ps,zs,.5,10,.3);
% 
% count = 1;
% tempp = [0,0,0];
% tempz = [0,0,0];
% tempd = [0,0,0];
% tempn = [0,0,0];
% jdelim = 0;
% for k = 1:length(ps)-1
%     if hs(k) > 0 && hs(k+1) < 0 && count == 1
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),eps(j),4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%         count = count+1;
%     elseif hs(k) < 0 && hs(k+1) > 0 && count == 1
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),eps(j),4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%         count = count+1;
%     elseif hs(k) > 0 && hs(k+1) < 0 && count == 2
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),eps(j),4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%         count = count+1;
%     elseif hs(k) < 0 && hs(k+1) > 0 && count == 2
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),eps(j),4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%         count = count+1;
%     elseif hs(k) > 0 && hs(k+1) < 0 && count == 3
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),eps(j),4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%     elseif hs(k) < 0 && hs(k+1) > 0 && count == 3
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),eps(j),4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,.3);
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%     end
% end
% %[temppN Ip] = sort(-tempp);
% %[tempp, I] = sort(tempp);
% %tempp = tempp;
% 
% tempp(tempp == 0) = NaN;
% pp(:,j) = tempp;
% 
% %tempz = tempz(I);
% tempz(tempz == 0) = NaN;
% zz(:,j) = tempz;
% 
% %tempd = tempd(I);
% tempd(tempd == 0) = NaN;
% dd(:,j) = tempd;
% 
% %tempn = tempn(I);
% tempn(tempn == 0) = NaN;
% nn(:,j) = tempn;
% 
% end
% 
% 
% % pp(:,1:idp(1)-1) = flipud(pp(:,1:idp(1)-1));
% % pp(:,idp(end)+1:end) = flipud(pp(:,idp(end)+1:end));
% % zz(:,idp(end)+1:end) = flipud(zz(:,idp(end)+1:end));
% % nn(:,idp(end)+1:end) = flipud(nn(:,idp(end)+1:end));
% % dd(:,idp(end)+1:end) = flipud(dd(:,idp(end)+1:end));
% ii = find(pp(1,:) > 0);
% if length(ii) > 0
%     for j = ii+1:length(eps)
%         temp1= pp(:,j);
%         [temp1 Ipp] = sort(temp1);
%         pp(:,j) = temp1;
%         temp2 = zz(:,j);
%         temp2 = temp2(Ipp);
%         zz(:,j) = temp2;
%         temp3 = dd(:,j);
%         temp3 = temp3(Ipp);
%         dd(:,j) = temp3;
%         temp4 = nn(:,j);
%         temp4 = temp4(Ipp);
%         nn(:,j) = temp4;
%         
%     end
% end
% 
% idxz1 = find(zz(1,:) <0);
% idxz3 = find(zz(3,:) < 0);
% 
% ii2 = find(pp(2,:)>0);
% % if length(ii2) > 0
% %    for j= ii2(end)+1:length(eps)
% %     pp(:,j) = flipud(pp(:,j));
% %     zz(:,j) = flipud(zz(:,j));
% %     nn(:,j) = flipud(nn(:,j));
% %     dd(:,j) = flipud(dd(:,j));
% %    end
% % end
% idp = find(pp(2,:) > 0);
% pp(:,idp(end)+1:end) = flipud(pp(:,idp(end)+1:end));
% nn(:,idp(end)+1:end) = flipud(nn(:,idp(end)+1:end));
% zz(:,idp(end)+1:end) = flipud(zz(:,idp(end)+1:end));
% dd(:,idp(end)+1:end) = flipud(dd(:,idp(end)+1:end));
% figure
% hold on
% plot(eps,pp,'b','linewidth',2)
% plot(eps,zz,'r','linewidth',2)
% plot(eps,dd,'k','linewidth',2)
% plot(eps,nn,'Color',color,'linewidth',2)
% 
% if length(ii2) > 0
%     xx = [eps(ii2(1)-1),eps(ii2(end)+1)];
%     x = [xx,fliplr(xx)];
%     ylb = [0,0];
%     yub = [14,14];
%     inBetween = [ylb,fliplr(yub)];
%     hh = fill(x, inBetween,'k','LineStyle','-','Marker','none');
%     set(hh,'facealpha',.3)
%     inBetween = [ylb,fliplr(yub)];
% end
% if length(idxz1) > 0
% xlim([0 eps(idxz1(1)-1)])
% end
% if length(idxz3)>0
% xlim([eps(idxz3(end)+1) eps(end)])
% end
% ylim([0 inf])
% xlabel('$\tilde{\alpha}$','fontsize',18)
% hold off
% 
% 
% %5%%%%%%%%%%%%555555555555555555555555555
% 
% clear all
% set(0,'defaultTextInterpreter','latex');
% ps = linspace(.000001,10,50000);
% 
% zstar = @(p,alpha,xi,a) ((p.^2./(1+p.^2))-((alpha.*p)./(xi+p))).*(1/a);
% dstar = @(p,z,gamma,psi,epsilon) (1./psi).*((1-gamma).*(p.^2./(1+p.^2)).*z + epsilon.*p);
% n_b = @(p,d,c,theta,k,psi,s) (p.*(1-p./c) + s.*(k-theta)-psi.*d);
% nstar = @(p,d,c,theta,k,psi,s) (-n_b(p,d,c,theta,k,psi,s) + sqrt((n_b(p,d,c,theta,k,psi,s)).^2+4.*s.*k.*(s.*theta+psi.*d)))*(1/(2.*s));
% 
% f = @(n,p,k,c) (n./(k+n)).*p.*(1-p./c);
% g = @(p,z,epsilon) (p.^2./(1+p.^2)).*z + epsilon.*p;
% h = @(n,p,z,k,c,epsilon) f(n,p,k,c) - g(p,z,epsilon);
% 
% 
% % toxicity/hyposixa half saturation 
% eps = linspace(.001,2,1000);
% for j = 1:length(eps)
% zs = zstar(ps,1.35,4.5,.4);
% ds = dstar(ps,zs,.5,.15,eps(j));
% ns = nstar(ps,ds,10,4,.5,.15,.3);
% hs = h(ns,ps,zs,.5,10,eps(j));
% 
% count = 1;
% tempp = [0,0,0];
% tempz = [0,0,0];
% tempd = [0,0,0];
% tempn = [0,0,0];
% jdelim = 0;
% for k = 1:length(ps)-1
%     if hs(k) > 0 && hs(k+1) < 0 && count == 1
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),1.35,4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,eps(j));
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%         count = count+1;
%     elseif hs(k) < 0 && hs(k+1) > 0 && count == 1
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),1.35,4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,eps(j));
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%         count = count+1;
%     elseif hs(k) > 0 && hs(k+1) < 0 && count == 2
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),1.35,4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,eps(j));
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%         count = count+1;
%     elseif hs(k) < 0 && hs(k+1) > 0 && count == 2
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),1.35,4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,eps(j));
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%         count = count+1;
%     elseif hs(k) > 0 && hs(k+1) < 0 && count == 3
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),1.35,4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,eps(j));
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%     elseif hs(k) < 0 && hs(k+1) > 0 && count == 3
%         tempp(count) = ps(k);
%         tempz(count) = zstar(ps(k),1.35,4.5,.4);
%         tempd(count) = dstar(ps(k),tempz(count),.5,.15,eps(j));
%         tempn(count) = nstar(ps(k),tempd(count),10,4,.5,.15,.3);
%     end
% end
% [tempp, I] = sort(tempp);
% tempp(tempp == 0) = NaN;
% pp(:,j) = tempp;
% 
% tempz = tempz(I);
% tempz(tempz == 0) = NaN;
% zz(:,j) = tempz;
% 
% tempd = tempd(I);
% tempd(tempd == 0) = NaN;
% dd(:,j) = tempd;
% 
% tempn = tempn(I);
% tempn(tempn == 0) = NaN;
% nn(:,j) = tempn;
% 
% end
% ii = find(pp(1,:) > 0);
% if length(ii) > 0
%     for j = ii+1:length(eps)
%         temp1= pp(:,j);
%         [temp1 Ipp] = sort(temp1);
%         pp(:,j) = temp1;
%         temp2 = zz(:,j);
%         temp2 = temp2(Ipp);
%         zz(:,j) = temp2;
%         temp3 = dd(:,j);
%         temp3 = temp3(Ipp);
%         dd(:,j) = temp3;
%         temp4 = nn(:,j);
%         temp4 = temp4(Ipp);
%         nn(:,j) = temp4;
%         
%     end
% end
% 
% idxz1 = find(zz(1,:) <0);
% idxz3 = find(zz(3,:) < 0);
% 
%  ii2 = find(pp(2,:)>0);
% % if length(ii2) > 0
% %    for j= ii2(end)+1:length(eps)
% %     pp(:,j) = flipud(pp(:,j));
% %     zz(:,j) = flipud(zz(:,j));
% %     nn(:,j) = flipud(nn(:,j));
% %     dd(:,j) = flipud(dd(:,j));
% %    end
% % end
% figure
% hold on
% plot(eps,pp(1,:),'b','linewidth',2)
% plot(eps,pp(2,:),'b','linewidth',2)
% plot(eps,pp(3,:),'b','linewidth',2)
% plot(eps,zz,'r','linewidth',2)
% plot(eps,dd,'k','linewidth',2)
% plot(eps,nn,'Color',color,'linewidth',2)
% 
% if length(ii2) > 0
%     xx = [eps(ii2(1)-1),eps(ii2(end)+1)];3
%     x = [xx,fliplr(xx)];
%     ylb = [0,0];
%     yub = [14,14];
%     inBetween = [ylb,fliplr(yub)];
%     hh = fill(x, inBetween,'k','LineStyle','-','Marker','none');
%     set(hh,'facealpha',.3)
%     inBetween = [ylb,fliplr(yub)];
% end
% if length(idxz1) > 0
% xlim([0 eps(idxz1(1)-1)])
% end
% if length(idxz3)>0
% xlim([eps(idxz3(end)+1) eps(end)])
% end
% ylim([0 inf])
% xlabel('$\tilde{\epsilon}$','fontsize',18)
% hold off
% 
% 
% 
% 
% 
% 
xis = linspace(.05,7.5,500);
alphas = linspace(0.05,3,500);
[X A] = meshgrid(xis, alphas);
for j = 1:500
    for jj = 1:500
    R0(j,jj) = (10.*(1-1/R0p)*(X(j,jj)+10.*(1-1/R0p)))./(A(j,jj).*(1+10.^2.*(1-1/R0p)^2));
    end
end
levels = 1:2*((max(max(R0))-1)./18):max(max(R0))
cmap = colormap(autumn(4));
cmap = [[0,0,1];cmap];
figure
hold on
colormap(flipud(autumn))
contourf(X,A,R0,[1,levels(2:end-1)])
hold off
set(gca,'FontSize',18)


end


