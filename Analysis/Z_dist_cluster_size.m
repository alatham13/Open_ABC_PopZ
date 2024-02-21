clear all;
close all;
step=3000; % number of timesteps (after removing equilibrium steps)

% Divide the length of the box (in nm) by the number of bins
dZ=500/100;

% protein density by Z-positon -------------------------------------------------------------
FL150=importdata('protein_hist_mindist.txt',' ',1);
FL150=FL150.data;

% Normalize
FL150(:,2)=FL150(:,2)./(step);
% Normalize the mass per step by the volume of each bin (in mL)
% xy is the xy box length taken from npt.xml
xy=29.263652792064985;
dV=dZ*xy*xy/(10^(21));
% Convert dalton to mg
FL150(:,2)=FL150(:,2)*1.6605*10^(-21);
FL150(:,2)=FL150(:,2)./(dV);

figure;
hold on;
plot(FL150(:,1)/10,FL150(:,2),'Linewidth',6);
set(gca,'FontSize',52,'FontName','Helvetica','Linewidth',4);
legend({'FL-PopZ 150 mM'},'location','northeast','FontSize',52,'FontName','Helvetica');
axis([-250 250 0 350]);
box on;


% cluster size -------------------------------------------------------------
FL=importdata('cluster_size_mindist.txt');

% take histogram. Maximum value is 100 trimers per cluster
edges=[-0.5:1:100.5];
[N1,edges]=histcounts(FL(:,1),edges,'Normalization','pdf');
edges2=[0:1:100];

figure;
hold on;
plot(edges2,N1,'Linewidth',6);
set(gca,'FontSize',52,'FontName','Helvetica','Linewidth',4);
legend({'FL-PopZ 150 mM'},'location','northeast','FontSize',52,'FontName','Helvetica');
axis([0 100 0 0.5]);
box on;
