function []= getCont
n0s = linspace(.1,14,300);
for j = 1:300
    fprintf(' \n')
    fprintf('ITERATION %i of %i\n',j,300)
    fprintf(' \n')
    EquibsC5(j,:) = FindEquibs(n0s(j),.5);
end
save('EquibsC5.mat','EquibsC5')
end