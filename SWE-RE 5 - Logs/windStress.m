%% WindStress

clc
clear all
close all


rhoA=1.1;

u=0:0.5:60;

cd=(-0.016*u.*u + 0.967*u + 8.058)*0.0001;

for i=1:size(cd,2)
    if(u(i)>50.33)
        cd(i)=4.1047/(u(i)^2);
    end
end

tau=cd.*u.*u*rhoA;

cdB=(0.067*u + 0.75)*0.001;
for i=1:size(cdB,2)
    if(cdB(i)>0.003)
        cdB(i)=.003;
    end
end

tauB=cdB.*u.*u*rhoA;

figure(1)
subplot(2,1,1)
hold on
plot(u,cd)
plot(u,cdB)
legend('SWE','ADCIRC')
title('Drag Coeff')
xlabel('Wind Speed m/s')
subplot(2,1,2)
hold on
plot(u,tau)
plot(u,tauB)
legend('SWE','ADCIRC')
title('WindStress')
xlabel('Wind Speed m/s')
