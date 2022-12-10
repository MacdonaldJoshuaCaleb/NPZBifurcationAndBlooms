function [] = MakeFigures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make manuscript figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

its = [.5,.85,1.6,3.2];
s = strings([length(its),1]);
s(1) = 'Favor zoopkankton';
s(2) = 'Ideal';
s(3) = 'Harmful algal blooms';
s(4) = 'Zooplankton extinction';
for j = 1:4
fig = ScenarioPlots(its(j),s(j));
str = string(strcat('Fig_',num2str(j)))
baseFileName = sprintf(str);
fname = '/home/anon/Documents/MATLAB/NPZD/FiguresForManuscript';
saveas(gca, fullfile(fname, baseFileName), 'png');
close all
end

end