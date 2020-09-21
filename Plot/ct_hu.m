% This script plots data from the ct_hu.cpp test program
clear;

% Always assume the data is in the "Build" folder on the same level as the
% "Plot" directory.
file_name = '../Build/ct_hu.txt';

% Read in data
fid = fopen(file_name);
C = textscan(fid, "%s %f %f %d %d","delimiter",",");
fclose(fid);

names = C{1};
density = C{2};
zeff = C{3};
hu_ave = (C{4} + C{5}) ./ 2;

% Sort based on density and HU
[sort_density, density_ind] = sort(density);
[sort_zeff, zeff_ind] = sort(zeff);
[sort_hu, hu_ind] = sort(hu_ave);
for n = 1:length(names)
  density_names{n} = names{density_ind(n)};
  zeff_names{n} = names{zeff_ind(n)};
  hu_names{n} = names{hu_ind(n)};
end

% Display bar graphs of density and HU
back_color = [0.15 0.15 0.15];

figure('position', [100 100 1000 800])
bar(sort_density, 'c')
ylabel('Density (g/cm^{3})')
set(gca, 'xticklabel', density_names)
set(gca, "fontsize", 20);
set(gca, 'color', back_color)
set(gca, 'xcolor', 'w')
set(gca, 'ycolor', 'w')
set(gcf, 'color', back_color)

figure('position', [120 120 1000 800])
bar(sort_zeff, 'c')
ylabel('Z_{eff}')
set(gca, 'xticklabel', zeff_names)
set(gca, "fontsize", 20);
set(gca, 'color', back_color)
set(gca, 'xcolor', 'w')
set(gca, 'ycolor', 'w')
set(gcf, 'color', back_color)

figure('position', [140 140 1000 800])
bar(sort_hu, 'c')
ylabel('HU')
set(gca, 'xticklabel', hu_names)
set(gca, "fontsize", 20);
set(gca, 'color', back_color)
set(gca, 'xcolor', 'w')
set(gca, 'ycolor', 'w')
set(gcf, 'color', back_color)