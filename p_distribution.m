clear all;
close all;
clc;

para = [30 0.05 18/5 3/5];

fnum=2;

para1=importdata(['ipara1_', num2str(fnum), '.mat']);
para2=importdata(['ipara2_', num2str(fnum), '.mat']);
para3=importdata(['ipara3_', num2str(fnum), '.mat']);
para4=importdata(['ipara4_', num2str(fnum), '.mat']);

a1=importdata(['apara1_', num2str(fnum), '.mat']);
a2=importdata(['apara2_', num2str(fnum), '.mat']);
a3=importdata(['apara3_', num2str(fnum), '.mat']);
a4=importdata(['apara4_', num2str(fnum), '.mat']);

sample=importdata(['sampletra_', num2str(fnum), '.mat']);
sample=[0,sample];

burn=200;
prod=para1(burn:end);
deg=para2(burn:end);
alpha=para3(burn:end);
beta=para4(burn:end);

timespan=100;
true_tra=mean_trajectory([0:timespan],[para(1) para(2) para(3) 1/para(4)]);
infer_tra=mean_trajectory([0:timespan],[mean(prod),mean(deg),mean(alpha),mean(beta)]);
% infer_tra=mean_trajectory([0:timespan],[prod(end),deg(end),alpha(end),beta(end)]);


subplot(2,2,1)
scatter(prod,deg); hold on;
scatter([para(1)],[para(2)],'filled'); hold off;
subplot(2,2,2)
scatter(alpha, 1./beta); hold on;
scatter([para(3)],[para(4)],'filled'); hold off;
subplot(2,2,3)
scatter(alpha.*beta,alpha.*(beta.^2)); hold on;
scatter([para(3)/para(4)],[para(3)/(para(4)^2)],'filled'); hold off;
subplot(2,2,4)
plot(0:100,sample,'o'); hold on;
plot(0:100,true_tra,'b-');
plot(0:100,infer_tra,'r-'); hold off;

% 
% inferred=[ipara1(end),ipara2(end),ipara3(end),1/ipara4(end)];
% infer_tra=mean_trajectory([0;timespan],[inferred(1),inferred(2),inferred(3),1/inferred(4)]);









