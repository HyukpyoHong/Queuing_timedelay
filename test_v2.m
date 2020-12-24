clear all;
close all;
clc;

fnum=2;

rng(fnum);

T_max = 100;
Iter_max = 1000000;

para = [35.4 0.015 5.89 0.89];
para = [30 0.05 18/5 3/5];

tprod=para(1);
tdeg=para(2);
talpha=para(3);
tbeta=para(4);

sampletra=[];
ipara1=[];
ipara2=[];
ipara3=[];
ipara4=[];
apara1=[];
apara2=[];
apara3=[];
apara4=[];

theta=[tprod tdeg talpha 1/tbeta];
theta_init = [(1/2)*tprod (1/2)*tdeg (1/2)*talpha 1/((1/2)*tbeta)];

[Xt, tspan, Xbirth, Xdeath] = Gillespie_delayX(theta, T_max, Iter_max);
data = cumsum(Xbirth - Xdeath);
timespan = (1:T_max)';

var_list = data + 1;
prior_mean = 1/2 * theta;
prior_var = [1000, 1, 10, 1];

burn = 0;
thin_rate = 1;
effnum = 250;

[rr, aa] = MCMC_function(data, timespan, var_list, prior_mean, prior_var, theta_init, burn, thin_rate, effnum);

%cc = clock;
%timestamp = [num2str(cc(1)) num2str(cc(2),'%02d') num2str(cc(3),'%02d') num2str(cc(4),'%02d') num2str(cc(5),'%02d') num2str(floor(cc(6)),'%02d')];
%save(['Queueing_MCMC_' timestamp]);

sampletra=[sampletra; data'];
ipara1=[ipara1; rr(:,1)'];
ipara2=[ipara2; rr(:,2)'];
ipara3=[ipara3; rr(:,3)'];
ipara4=[ipara4; rr(:,4)'];

apara1=[apara1; aa(:,1)'];
apara2=[apara2; aa(:,2)'];
apara3=[apara3; aa(:,3)'];
apara4=[apara4; aa(:,4)'];

save(['sampletra_', num2str(fnum), '.mat'],'sampletra');

save(['ipara1_', num2str(fnum), '.mat'],'ipara1');
save(['ipara2_', num2str(fnum), '.mat'],'ipara2');
save(['ipara3_', num2str(fnum), '.mat'],'ipara3');
save(['ipara4_', num2str(fnum), '.mat'],'ipara4');

save(['apara1_', num2str(fnum), '.mat'],'apara1');
save(['apara2_', num2str(fnum), '.mat'],'apara2');
save(['apara3_', num2str(fnum), '.mat'],'apara3');
save(['apara4_', num2str(fnum), '.mat'],'apara4');

%%

inferred=[ipara1(end),ipara2(end),ipara3(end),1/ipara4(end)];
infer_tra=mean_trajectory([0;timespan],[inferred(1),inferred(2),inferred(3),1/inferred(4)]);

plot(0:1:100,[0,sampletra],'-o'); hold on;
plot(0:1:100,infer_tra','-r'); hold off;






