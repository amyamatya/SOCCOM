%% Temperature vs. Cphyto Anomaly 
%0 = 6/2006, 165 = 3/2019 (length 166)
hold on
plot(1:166,[repmat(bbpttl55,1,13) bbpttl55(1:10)], 'color',[.87 .87 .87], 'linewidth',5);
plot(71:166, floatplot55,'c','linewidth',1);
plot(1:144, dayplot55,'--b','linewidth',2); 
plot(1:144, nightplot55,':r','LineWidth',2);
hold off
axis tight
legend('Overall','Float','Satellite (day)','Satellite (night)');
xlabel('date')
ylabel('bbp');
title('monthly bbp averages 60\circ (2006 - 2019)');
xticks([1:12:166]);
xlab = 2006:2019;
xticklabels(xlab);
axis tight;

%6/2006 to 4/2017 

%% Temp vs. Cphyto, filled
load sst_anomalies.mat;
firstidx = find(dates == datenum('01-Jun-2006'));
lastidx = find(dates == datenum('01-Apr-2017'));
sstdata = data(firstidx:lastidx);
diffnan = geodiv(11).geodiff(6:136);%(~isnan(daydiff55));

fig10 = figure(10);
hold on
yyaxis left
%line([1 137],[0 0],'Color',[.9 .9 .9],'linewidth',5,'Handlevisibility','off');  %length no nan = 131
a1 = area(daydiff55);a1.FaceColor = [.8 .4 .4]; a1.FaceAlpha = 0.5;
ylim([-.0006 .0014]);
ylabel('Monthly b_{bp} Anomaly, Averaged 2006-2019  (m^{-1})');

yyaxis right
a2 = area(5:135,sstdata); a2.FaceColor = [.4 .6 .8]; a2.FaceAlpha = 0.5;
legend('b_{bp} Anomaly, below 55\circS','Surface Temperature Anomaly below 55\circS','Location','NorthWest');
title('Temperature and b_{bp} Monthly Anomaly, 55\circS, (2006-2017)')
ylabel('Monthly Temperature Anomaly, Averaged 1982-2019 (C\circ)');
xlim([5 137]);
xticks([1:12:166]); xlab = 2006:2019; xticklabels(xlab);
xlabel('Year');
hold off
print(fig10, '/Users/aamatya/Desktop/Summer 2019/10-sstanomaly','-dpdf');
%% Temp vs. Cphyto, inverted
diffnan2 = diffnan; sstdata2 = sstdata;
nanidx = (isnan(diffnan2) == 1);
diffnan2(nanidx) = [];
sstdata2(nanidx) = [];
corrcoef(diffnan2,sstdata2)

fig11 = figure(11);
hold on
plot([4 137], [0 0], 'linewidth',5,'color',[.9 .9 .9]); alpha(.3w);
yyaxis left
plot(geodiv(11).geodiffphyto); % 6-137 are valued
ylim([-4 4]);
ylabel('C_{phyto} (m^{-1}), Averaged 2006-2019');
yyaxis right
plot(6:136,sstdata);
ylim([-.68 .68]);
ylabel('SST (C\circ), Averaged 1982-2019');
title('Temperature and C_{phyto} Monthly Anomaly Below 55\circS,  r^2 = -0.1333')
text(92,.63,'*Temperature inverted','Color',[.8 .45 .35],'fontsize',13);
xlim([5 137]);
xticks([1:12:166]); xlab = 2006:2019; xticklabels(xlab);
xlabel('Year (2006-2017)');
axis ij
hold off
print(fig11, '/Users/aamatya/Desktop/Summer 2019/11-anomalytrs','-dpdf');

%%
firstidx = find(dates == datenum('01-Jun-2006'));
lastidx = find(dates == datenum('01-Apr-2017'));
sstdata = data(firstidx:lastidx);
diffnan = geodiv(11).geodiff(6:136);

hold on
yyaxis right
plot(6:136,sstdata, 'color',[.8 .4 .2]);
temp = sstdata;
id = find(temp >= 0); temp(id) = 0;
a2 = area(6:136, temp);
 a2.FaceColor = [.8 .4 .2]; a2.FaceAlpha = 0.4;
ylim([-.68 .68]);
ylabel('SST (C\circ), Averaged 1982-2019');
title('Temperature and b_{bp} Monthly Anomaly Below 55\circS,  r^2 = -0.1333')
text(92,.63,'*Temperature inverted','Color',[.8 .45 .35],'fontsize',13);
xlim([5 137]);
xticks([1:12:166]); xlab = 2006:2019; xticklabels(xlab);
xlabel('Year (2006-2017)');
axis ij

% plot([5 137], [0 0], 'linewidth',5,'color',[.9 .9 .9]);
yyaxis left
temp = geodiv(11).geodiff;
plot(geodiv(11).geodiff, 'color',[.5 .5 .5]); % 6-137 are valued
id2 = find(temp < 0);
temp(id2) = 0;
a1 = area(temp);
a1.FaceColor = [.5 .5 .5]; a1.FaceAlpha = 0.4;
ylim([-.0014 .0014]);
ylabel('b_{bp} (m^{-1}), Averaged 2006-2019');
hold off

figure,
yyaxis right
a2 = area(6:136,sstdata);
a2.FaceColor = [.8 .4 .2]; a2.FaceAlpha = 0.4;
ylim([-.68 .68]);
ylabel('SST (C\circ), Averaged 1982-2019');
title('Temperature and b_{bp} Monthly Anomaly Below 55\circS,  r^2 = -0.1333')
text(92,.63,'*Temperature inverted','Color',[.8 .45 .35],'fontsize',13);
xlim([5 137]);
xticks([1:12:166]); xlab = 2006:2019; xticklabels(xlab);
xlabel('Year (2006-2017)');
axis ij

hold on
yyaxis left
a1 = area(geodiv(11).geodiff); % 6-137 are valued
a1.FaceColor = [.2 .2 .2]; a1.FaceAlpha = 0.4;
ylim([-.0014 .0014]);
ylabel('b_{bp} (m^{-1}), Averaged 2006-2019');
% plot([5 137], [0 0], 'linewidth',5,'color',[.9 .9 .9]);
hold off



