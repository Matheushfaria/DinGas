% load('x.mat');
% load('Machf_bell.mat');
% load('press_bell.mat');
% load('dens_bell.mat');
% load('Temp_bell.mat');
% load('Veloc_bell.mat');
% 
% Mach_bell=Machf(2:end);
% pressure_bell=press(2:end);
% density_bell=dens(2:end);
% temperature_bell=Temp(2:end);
% velocity_bell=Veloc(2:end);

% anynan(density_bell)
% anynan(Mach_bell)
% anynan(pressure_bell)
% anynan(temperature_bell)
 anynan(velocity_bell)

close all;
figure(1)
plot(x, Veloc(2:end))
figure(2)
plot(x, velocity_bell)



% k=1668-1595;
% 
% diff=velocity_bell(1668)-velocity_bell(1595);
% 
% acr=diff/k;
% 
% for i=1596:1667
%     velocity_bell(i)=velocity_bell(i-1)+acr;
% end

