function ret = NewEquibs(n0,c,ghat,s)
%close all
%n0 = 12;   % intial biomass
%s = .3;    % nutrient loss rate
%c = 10;    % phytoplankton carrying capacity
g = .29;   % zooplankton assimilation 
k = .5;    % half saturation nutrient uptake 
%ghat = .5; %.8 % assimilation scaled by ratio of phytoplankton growth and max zoo grazing
aa = linspace(.001,.999,250);
zs = zeros(3,length(aa));
ps = zeros(3,length(aa));
ns = zeros(3,length(aa));
pstar = @(z,a) sqrt((a.*z)./(1-a.*z));
nstar = @(z,a) (1/s).*(s.*n0 - (g.*a.*z).*z);
J = @(n,p,z,a) [(2*p^3*z)/(p^2 + 1)^2 - (n*p)/(c*k + c*n) + (n*(c - p))/(c*k + c*n) - (2*p*z)/(p^2 + 1), -p^2/(p^2 + 1), (p*(c - p))/(c*k + c*n) - (c*n*p*(c - p))/(c*k + c*n)^2;
              (2*ghat*p*z)/(p^2 + 1) - (2*ghat*p^3*z)/(p^2 + 1)^2, (ghat*p^2)/(p^2 + 1) - 2*a*ghat*z,0;
              (n*p)/(c*k + c*n) - (n*(c - p))/(c*k + c*n) - (2*p*z*(g - 1))/(p^2 + 1) + (2*p^3*z*(g - 1))/(p^2 + 1)^2,-(p^2*(g - 1))/(p^2 + 1),(c*n*p*(c - p))/(c*k + c*n)^2 - (p*(c - p))/(c*k + c*n) - s];
for w = 1:length(aa)
a = aa(w); % zooplankton loss rate 
zz = linspace(0,min(1/a,sqrt((n0*s/g)*(1/a))),250);
zstar1 = @(z) (nstar(z,a)./(k+nstar(z,a))).*pstar(z,a);
zstar2 = @(z) (nstar(z,a)./(k+nstar(z,a))).*(pstar(z,a).^2./c)+a.*z.^2;
zstar = @(z) zstar1(z) - zstar2(z);
% figure
% plot(zz,(zstar1(zz)).^2,zz,(zstar2(zz)).^2,'linewidth',2)
% 
% figure
% plot(zz,zstar(zz),zz,zeros(length(zz),1),'k','linewidth',2)
% ylim([-1,1])
% xlim([0,zz(end)])

% for j = 1:3
%     for k = 1:length(zs)
%         if zs(j,k) == 0
%             zs(j,k) = NaN;
%             ps(j,k) = NaN;
%             ns(j,k) = NaN;
%         end
%     end
% end

root = zeros(3,1);
count = 0;
for j = 2:length(zz)
    if zstar(zz(j))>0 && zstar(zz(j-1))<0 && count == 0 
        root(1) = fzero(zstar,zz(j));
        if root(1) == NaN
            root(1) = (zstar(zz(j))+zstar(zz(j-1)))./2;
        end
        count = count+1;
        temp = j;
    end
    if zstar(zz(j))>0 && zstar(zz(j-1))<0 && count == 1 && temp < j
        root(2) = fzero(zstar,zz(j));
        if root(2) == NaN
            root(2) = (zstar(zz(j))+zstar(zz(j-1)))./2;
        end
        count = count+1;
        temp = j;
    end
    if zstar(zz(j))>0 && zstar(zz(j-1))<0 && count == 2 && temp < j
        root(3) = fzero(zstar,zz(j));
        if root(3) == NaN
            root(3) = (zstar(zz(j))+zstar(zz(j-1)))./2;
        end
        count = count+1;
        temp = j;
    end
    if zstar(zz(j))<0 && zstar(zz(j-1))>0 && count == 0 
        root(1) = fzero(zstar,zz(j));
        if root(1) == NaN
            root(1) = (zstar(zz(j))+zstar(zz(j-1)))./2;
        end
        count = count+1;
        temp = j;
    end
     if zstar(zz(j))<0 && zstar(zz(j-1))>0 && count == 1 && temp < j
        root(2) = fzero(zstar,zz(j));
        if root(2) == NaN
            root(2) = (zstar(zz(j))+zstar(zz(j-1)))./2;
        end
        count = count+1;
        temp = j;
     end
     if zstar(zz(j))<0 && zstar(zz(j-1))>0 && count == 2 && temp < j
        root(3) = fzero(zstar,zz(j));
        if root(3) == NaN
            root(3) = (zstar(zz(j))+zstar(zz(j-1)))./2;
        end
        count = count+1;
        temp = j;
        
    end
end
%root = root(root>0);
zs(:,w) = root;
%ps(:,w) = pstar(zs(:,w));
%ns(:,w) = nstar(zs(:,w));
%for w = 1:length(aa)
%    for j = 1:3
%        if ns(j,w) == n0
%            ns(j,w) = 0;
%        end
%    end
%end
% %root = root(root>0);
 if w > 1
    if min(zs(:,w-1))> 0 && min(zs(:,w))==0
        bound = w;
        %zs(:,bound:end) = sort(zs(:,bound:end)); 
    end
    if min(zs(:,w-1)) == 0 && length(root(root > 0)) == 2
        bound = w;
        %zs(:,bound:end) = sort(zs(:,bound:end)); 
    end
       
 end
end
bound
zs(:,bound:end) = sort(zs(:,bound:end)); 
%zs(:,bound:end) = sort(zs(:,bound:end));
ps = pstar(zs,aa);
ns = nstar(zs,aa);
% ps(3,bound:end) = ps(1,bound:end);
% ps(1,bound:end) = zeros(1,length(ps(1,bound:end)));
% % zs(3,lbound) = (zs(2,lbound)+zs(3,lbound+1))./2;
% % ps(3,lbound) = (ps(2,lbound)+ps(3,lbound+1))./2;
% zs(1,bound) = (zs(1,bound-1)+zs(2,bound-1))./2;
% zs(2,bound) = (zs(1,bound-1)+zs(2,bound-1))./2;
% ps(1,bound) = (ps(1,bound-1)+ps(2,bound-1))./2;
% ps(2,bound) = (ps(1,bound-1)+ps(2,bound-1))./2;
% zs(3,lbound-1) = (zs(2,lbound)+zs(3,lbound))./2;
% zs(2,lbound-1) = (zs(2,lbound)+zs(3,lbound))./2;
% ps(3,lbound-1) = (ps(2,lbound)+ps(3,lbound))./2;
% ps(2,lbound-1) = (ps(2,lbound)+ps(3,lbound))./2;
 for w = 1:length(ns(1,:))
      for j = 1:3
         if ns(j,w) == n0
             ns(j,w) = 0;
         end
     end
 end
% ns(3,bound:end) = ns(1,bound:end);
% ns(1,bound:end) = zeros(1,length(ns(1,bound:end))); 
% ns(3,lbound) = (ns(2,lbound)+ns(3,lbound+1))./2;
% ns(1,bound) = (ns(1,bound-1)+ns(2,bound-1))./2;
% ns(2,bound) = (ns(1,bound-1)+ns(2,bound-1))./2;
% ns(3,lbound-1) = (ns(2,lbound)+ns(3,lbound))./2;
% ns(2,lbound-1) = (ns(2,lbound)+ns(3,lbound))./2;
 for w = 1:length(zs(1,:))
     for j = 1:3
        if zs(j,w) == 0
            zs(j,w) = NaN;
            ps(j,w) = NaN;
             ns(j,w) = NaN;
         end
     end
 end
 
temp = zs(2,:);
if length(temp(temp>0)) == 1
    zs(1,bound) = NaN;
    zs(2,bound) = NaN;
    ps(1,bound) = NaN;
    ps(2,bound) = NaN;
    ns(1,bound) = NaN;
    ns(2,bound) = NaN;
    zs(1,bound:end) = zs(3,bound:end);
    zs(3,bound:end) = zs(2,bound:end);
    ps(1,bound:end) = ps(3,bound:end);
    ps(3,bound:end) = ps(2,bound:end);
    ns(1,bound:end) = ns(3,bound:end); 
    ns(3,bound:end) = ns(2,bound:end);
end

H1 = zeros(3,length(aa));
H2 = zeros(3,length(aa));
H3 = zeros(3,length(aa));
psstab = zeros(3,length(aa));
zsstab = zeros(3,length(aa));
nsstab = zeros(3,length(aa));
DH2 = zeros(3,length(aa));

for w = 1:length(zs(1,:))
     for j = 1:3
        H1(j,w) = NaN;
        H2(j,w) = NaN;
        H3(j,w) = NaN;
        psstab(j,w) = NaN;
        zsstab(j,w) = NaN;
        nsstab(j,w) = NaN;
        DH2(j,w) = NaN;
        if zs(j,w) > 0
        M = J(ns(j,w),ps(j,w),zs(j,w),aa(w));
        cofs = charpoly(M);
        H1(j,w) = cofs(2);
        H2(j,w) = cofs(2)*cofs(3) - cofs(4);
        H3(j,w) = cofs(4)*H2(j,w);
        if H1(j,w) > 0 && H2(j,w) > 0 && H3(j,w) > 0
            zsstab(j,w) = zs(j,w);
            psstab(j,w) = ps(j,w);
            nsstab(j,w) = ns(j,w);
        end
            if w > 1 && H2(j,w) ~= 0 && H2(j,w-1) ~= 0
                DH2(j,w) = (H2(j,w) - H2(j,w-1))/(aa(w)-aa(w-1));
            end
        end
     end
end
index1 = 0;
index2 = 0;
its = length(aa);
for k = 2:length(aa)
    if H1(1,k) > 0
        if H2(1,k-1) < 0 && H2(1,k) > 0 && DH2(1,k) ~= 0
                index1 = k;
        end
        if H2(1,k-1) > 0 && H2(1,k) < 0 && DH2(1,k) ~= 0
                index2 = k;
        end
    end
end
if index1 > 0 && index2 > 0
lb1 = min(index1,index2)
ub1 = max(index1,index2)
end

index3 = 0;
index4 = 0;

for k = 2:length(aa)
    if H1(2,k) > 0
        if H2(2,k-1) < 0 && H2(2,k) > 0 && DH2(2,k) ~= 0
                index3 = k;
        end
        if H2(2,k-1) > 0 && H2(2,k) < 0 && DH2(2,k) ~= 0
                index4 = k;
        end
    end
end
if index3 > 0 && index4 > 0
lb2 = min(index3,index4)
ub2 = max(index3,index4)
end

index5 = 0;
index6 = 0;

for k = 2:its
    if H1(3,k) > 0
        if H2(3,k-1) < 0 && H2(3,k) > 0 && DH2(3,k) ~= 0
                index5 = k;
        end
        if H2(3,k-1) > 0 && H2(3,k) < 0 && DH2(3,k) ~= 0
                index6 = k;
        end
    end
end
if index5 > 0 && index6 > 0
lb3 = min(index5,index6)
ub3 = max(index5,index6)
end



figure
hold on
plot(aa,psstab,'c',aa,ps,'c--','linewidth',2)
plot(aa,zsstab,'r',aa,zs,'r--','linewidth',2)
plot(aa,nsstab,'y',aa,ns,'y--','linewidth',2)
if index1 > 0 && index2 > 0
plot(aa(lb1:ub1),ps(1,lb1:ub1),'c*','linewidth',2)
plot(aa(lb1:ub1),zs(1,lb1:ub1),'r*','linewidth',2)
plot(aa(lb1:ub1),ns(1,lb1:ub1),'y*','linewidth',2)
end
if index3 > 0 && index4 > 0
plot(aa(lb2:ub2),ps(2,lb2:ub2),'c*','linewidth',2)
plot(aa(lb2:ub2),zs(2,lb2:ub2),'r*','linewidth',2)
plot(aa(lb2:ub2),ns(2,lb2:ub2),'y*','linewidth',2)
end
if index5 > 0 && index6 > 0
plot(aa(lb3:ub3),ps(3,lb1:ub3),'c*','linewidth',3)
plot(aa(lb3:ub3),zs(3,lb1:ub3),'r*','linewidth',3)
plot(aa(lb3:ub3),ns(3,lb1:ub3),'y*','linewidth',3)
end
plot(aa,c*ones(length(aa),1),'c--','linewidth',2)
plot(aa,zeros(length(aa),1),'c--','linewidth',2)
plot(aa,n0*ones(length(aa),1),'y--','linewidth',2)
ylabel('concentration','FontSize',16)
xlabel('a','FontSize',16)
hold off

ret = zeros(1,250);
for j = 1:250
if psstab(1,j) > 0 && psstab(3,j) > 0
    ret(j) = 3;
end
if psstab(1,j) > 0 && ret(j) ~= 3
    ret(j) = 1;
end
if psstab(3,j) > 0 && ret(j) ~= 3
    ret(j) = 2;
end

end

