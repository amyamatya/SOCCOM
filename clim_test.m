function clim_test % average all data for each month, average float/sat day/sat night per month per year, plot
%% 1 - load data
% Access directory
cd('/Users/aamatya/Desktop/SOCCOM');
load Datosint2.mat;

% Simplify time variables into single-rowed arrays
for i = 1:length(Datosint2)
    Datosint2(i).month_simple = Datosint2(i).month(1,:);
    Datosint2(i).year_simple = Datosint2(i).year(1,:);
end

% Extract desired variables
floatmonth = extractfield(Datosint2,'month_simple');
floatyear = extractfield(Datosint2, 'year_simple');
floatbbp = extractfield(Datosint2,'bbp_mld'); 
floatlon = extractfield(Datosint2,'lon');
floatqf = extractfield(Datosint2, 'lat_QF');
floatlat = extractfield(Datosint2, 'lat');

satmonth = extractfield(Datoslid_SO,'month_simple');
satbbp = extractfield(Datoslid_SO, 'bbp');
satmonthrep = month(extractfield(Datoslid_SO, 'time'));
satyearrep = year(extractfield(Datoslid_SO, 'time'));
satlat = extractfield(Datoslid_SO, 'lat');
satlon = extractfield(Datoslid_SO, 'lon');
%% 2 - calc averages

% Store variables by month
for n = 1:12
    a = find(floatmonth == n);
    clim(n).bbp = floatbbp(a);
    clim(n).lat = floatlat(a);
    clim(n).lon = floatlon(a);
    clim(n).qf = floatqf(a);
    clear a
end

% Combine float and satellite data into single array
for n = 1:12
    data_bbp = [];
    data_lat = [];
    data_lon = [];
    a = find(satmonth == n);
    for ii = 1:length(a)
        b = a(ii);
        pre_bbp = (Datoslid_SO(b).bbp);
        pre_lat = Datoslid_SO(b).lat;
        pre_lon = Datoslid_SO(b).lon;
        data_bbp = cat(1,data_bbp,pre_bbp);
        data_lat = cat(1,data_lat,pre_lat);
        data_lon = cat(1,data_lon,pre_lon);
    end
    clim(n).bbp_sat = data_bbp;
    clim(n).lat_sat = data_lat;
    clim(n).lon_sat = data_lon;
end

for n = 1:12
    clim(n).bbp_total = cat(1,clim(n).bbp_sat,clim(n).bbp');
    clim(n).lat_total = cat(1,clim(n).lat_sat,clim(n).lat');
    clim(n).lon_total = cat(1,clim(n).lon_sat,clim(n).lon');
end

clear i ii n a b pre_bbp pre_lat pre_lon data_lat data_lon data_bbp
%% 3 - Map plots float/satellite locations
filename = gunzip('gshhs_c.b.gz', tempdir);
shorelines = gshhs(filename{1});
delete(filename{1})
levels = [shorelines.Level];
land = (levels == 1);

fig1 = figure(1)
for n = 1:6
    display(n)
    sb(1) = subplot(4,3,n);
    hold on
    worldmap([-90 -20],[0 360])
    scatterm(clim(n).lat_sat,clim(n).lon_sat,0.5, [1 0 0],'filled'); %colormap(summer)
    scatterm(clim(n).lat,clim(n).lon,3,[0 0 1], 'filled'); %colormap(spring)
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
    title(n);
    hold off
    
end

set(fig1,'PaperUnits', 'centimeters', 'PaperSize', [50 25], 'PaperPosition', [0 0 50 25])
print(fig1, '/Users/aamatya/Desktop/Summer 2019/1-mapfig','-dpdf');

fig2 = figure(2)
for n = 7:12
    display(n)
    subplot(4,3,n-6);
    hold on
    worldmap([-90 -20],[0 360])
    scatterm(clim(n).lat_sat,clim(n).lon_sat,0.5, [1 0 0],'filled'); %colormap(summer)
    scatterm(clim(n).lat,clim(n).lon,3,[0 0 1], 'filled'); %colormap(spring)
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
    title(n);
    hold off
    
end
set(fig2,'PaperUnits', 'centimeters', 'PaperSize', [50 25], 'PaperPosition', [0 0 50 25])
print(fig2, '/Users/aamatya/Desktop/Summer 2019/2-mapfig','-dpdf');
%% 4 - monthly averages

% All data
for n = 1:12 
    bbpttl(n) = nanmean(clim(n).bbp_total);
end

% Satellite 2006-2017
for i = 1:12 %year
    for j = 1:12 %month
        thisdate = find((satmonthrep == j) & (satyearrep - 2005 == i));
        day = thisdate(1:2:end);
        night = thisdate(2:2:end);
        clim(j).satdaymean(i) = nanmean(satbbp(day));
        clim(j).satnightmean(i) = nanmean(satbbp(night));
    end
end

dayplot = [];
nightplot = [];
for i = 1:12 %10
    for j = 1:12
        thisday = clim(j).satdaymean(i);
        thisnight = clim(j).satnightmean(i);
        dayplot = [dayplot thisday];
        nightplot = [nightplot thisnight];
    end
end

clear thisdate thisday thisnight day night ;

% Float
for i = 1:8 
    for j = 1:12
        thisdate = find((floatmonth == j) & (floatyear - 2011 == i));
        clim(j).floatmean(i) = nanmean(floatbbp(thisdate));
    end
end

floatplot = [];
for i=1:8
    for j = 1:12
        thisday = clim(j).floatmean(i);
        floatplot = [floatplot thisday];
    end
end
clear i j n xlab thisday thisdate 
%% 5 - Monthly satellite bbp average
fig3 = figure(3); %0 = 6/2006, 165 = 3/2019 (length 166)
hold on
% plot(5:166,[bbpttl(6:12) repmat(bbpttl,1,13)], 
plot(1:166,[repmat(bbpttl,1,13) bbpttl(1:10)], 'color',[.87 .87 .87], 'linewidth',6);
% Uncomment the following line and legend line to add float to plot
%plot(71:166, floatplot,'c','linewidth',1); 
plot(1:144, dayplot,'--b','linewidth',2); 
plot(1:144, nightplot,':r','LineWidth',2);
%legend('Overall','Float','Satellite (day)','Satellite (night)');
legend('Overall','Satellite (day)','Satellite (night)'); xlim([5 138]);ylim([.0015 .0022]);
xlabel('Year'); ylabel('b_{bp}  (m^{-1})'); title('Monthly LIDAR b_{bp} Average (2006 - 2017)');
xticks([1:12:166]);
xlab = 2006:2017;
xticklabels(xlab);
hold off
print(fig3, '/Users/aamatya/Desktop/Summer 2019/3-monthly','-dpdf');
%% Monthly satellite bbp anomaly
floatdiff = [];
bbpavg = repmat(bbpttl,1,9);
for i = 1:length(floatplot)
    floatdiff = [floatdiff (floatplot(i)-bbpavg(i))];
end
daydiff = [];
nightdiff = [];
bbpavg = repmat(bbpttl, 1, 12);
for i = 1:length(dayplot)
    daydiff = [daydiff (dayplot(i) - bbpavg(i))];
    nightdiff = [nightdiff (nightplot(i) - bbpavg(i))];
end

fig4 = figure(4);
hold on
plot([1 166], [0 0],'color',[.9 .9 .9],'linewidth',5,'HandleVisibility','off');
% Uncomment the following line and legend line to add float to plot
%plot(71:166, floatdiff,'c','linewidth',1);
plot(1:144, daydiff,'--b','linewidth',2);
plot(1:144, nightdiff,':r','linewidth',2);
xticks([1:12:166]);
xlab = 2006:2019;
xticklabels(xlab);
%legend('Float','Satellite (day)','Satellite (night)');
legend('Satellite (day)', 'Satellite (night)'); ylim([-.0002 .0002]);xlim([5 137]);
xlabel('Year');
ylabel('Deviation from Monthly b_{bp} Average, 2006-2019  (m^{-1})');
title('Monthly LIDAR b_{bp} Anomaly (2006-2017)');
hold off
print(fig4, '/Users/aamatya/Desktop/Summer 2019/4-anomaly','-dpdf');

%% 7 - Cphyto monthly & anomaly
%cphyto = 0.19 * (3.12 * (10^4) * (bbp700) + 3) + 8.7;

floatphyto = 0.19 * (3.12 * (10^4) .* (floatplot) + 3) + 8.7;
dayphyto = 0.19 * (3.12 * (10^4) .* (dayplot) + 3) + 8.7;
nightphyto = 0.19 * (3.12 * (10^4) .* (nightplot) + 3) + 8.7;
totalphyto = 0.19 * (3.12 * (10^4) .* (bbpttl) + 3) + 8.7;

fig5 = figure(5); %0 = 6/2006, max = 2019
hold on
plot(1:166,[repmat(totalphyto,1,13) totalphyto(1:10)], 'color',[.87 .87 .87], 'linewidth',5);
plot(71:166, floatphyto,'c','linewidth',1);
plot(1:144, dayphyto,'--b','linewidth',2); 
plot(1:144, nightphyto,':r','LineWidth',2);
hold off
axis tight
legend('Overall','Float','Satellite (day)','Satellite (night)');
xlabel('Year')
ylabel('mg C m^{-3}');
title('Monthly C_{phyto} Averages (2006 - 2019)');
xticks([1:12:166]);
xlab = 2006:2019;
xticklabels(xlab);
axis tight;
print(fig5, '/Users/aamatya/Desktop/Summer 2019/5-monthlyphyto','-dpdf');

floatdiffphyto = [];
bbpavgphyto = repmat(totalphyto,1,9);
for i = 1:length(floatphyto)
    floatdiffphyto = [floatdiffphyto (floatphyto(i)-bbpavgphyto(i))];
end

daydiffphyto = [];
nightdiffphyto = [];
bbpavgphyto = repmat(totalphyto, 1, 12);
for i = 1:length(dayphyto)
    daydiffphyto = [daydiffphyto (dayphyto(i) - bbpavgphyto(i))];
    nightdiffphyto = [nightdiffphyto (nightphyto(i) - bbpavgphyto(i))];
end

fig6 = figure(6);
hold on
plot([1 166], [0 0],'color',[.85 .85 .85],'linewidth',2,'HandleVisibility','off');
plot(71:166, floatdiffphyto,'c','linewidth',1);
plot(daydiffphyto,'--b','linewidth',2);
plot(nightdiffphyto,':r','linewidth',2);
axis tight
hold off
xticks([1:12:166]);
xlab = 2006:2019;
xticklabels(xlab);
legend('Float','Satellite (day)','Satellite (night)');
xlabel('Year');
ylabel('Deviation from Monthly Mean (mg C m^{-3})');
title('C_{phyto} Anomaly (2006-2019)');
print(fig6, '/Users/aamatya/Desktop/Summer 2019/6-anomalyphyto','-dpdf');

clear floatdiffphyto bbpavgphyto daydiffphyto nightdiffphyto dayphyto nightphyto floatphyto... 
    totalphyto xlab bbpavg daydiff nightdiff
%% 8  Limited to -60S
floatidx = find(floatlat <= -60);
satidx = find(satlat <= -60);

floatmonth60 = floatmonth(floatidx);
floatyear60 = floatyear(floatidx);
floatbbp60 = floatbbp(floatidx);

satmonthrep60 = satmonthrep(satidx);
satyearrep60 = satyearrep(satidx);
satbbp60 = satbbp(satidx);

for n = 1:12
    a = find(floatmonth60 == n);
    clim(n).d60bbp = floatbbp60(a);
    clear a
end

for n = 1:12
    b = find(satmonthrep60 == n);
    clim(n).d60bbp_sat = satbbp60(b);
end

for n = 1:12
    clim(n).d60bbp_total = cat(2, clim(n).d60bbp_sat, clim(n).d60bbp);
end

clear i ii n a b pre_bbp pre_lat pre_lon data_lat data_lon data_bbp

% All data
for n = 1:12 
    bbpttl60(n) = nanmean(clim(n).d60bbp_total);
end

% Satellite 2006-2017
for i = 1:12 %year
    for j = 1:12 %month
        thisdate = find((satmonthrep60 == j) & (satyearrep60 - 2005 == i));
        day = thisdate(1:2:end);
        night = thisdate(2:2:end);
        clim(j).d60satdaymean(i) = nanmean(satbbp60(day));
        clim(j).d60satnightmean(i) = nanmean(satbbp60(night));
    end
end

dayplot60 = [];
nightplot60 = [];
for i = 1:12 %10
    for j = 1:12
        thisday = clim(j).d60satdaymean(i);
        thisnight = clim(j).d60satnightmean(i);
        dayplot60 = [dayplot60 thisday];
        nightplot60 = [nightplot60 thisnight];
    end
end

clear thisdate thisday thisnight day night ;

% Float
for i = 1:8 
    for j = 1:12
        thisdate = find((floatmonth60 == j) & (floatyear60 - 2011 == i));
        clim(j).d60floatmean(i) = nanmean(floatbbp60(thisdate));
    end
end

floatplot60 = [];
for i=1:8
    for j = 1:12
        thisday = clim(j).d60floatmean(i);
        floatplot60 = [floatplot60 thisday];
    end
end
clear i j n xlab thisday thisdate 
%% Monthly bbp averages
fig7 = figure(7); %0 = 6/2006, 165 = 3/2019 (length 166)
hold on
plot(1:166,[repmat(bbpttl60,1,13) bbpttl60(1:10)], 'color',[.87 .87 .87], 'linewidth',5);
plot(71:166, floatplot60,'c','linewidth',1);
plot(1:144, dayplot60,'--b','linewidth',2); 
plot(1:144, nightplot60,':r','LineWidth',2);
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
print(fig7, '/Users/aamatya/Desktop/Summer 2019/7-monthly60','-dpdf');

%% Monthly bbp anomalies
floatdiff60 = [];
bbpavg60 = repmat(bbpttl60,1,9);
for i = 1:length(floatplot60)
    floatdiff60 = [floatdiff60 (floatplot60(i)-bbpavg60(i))];
end
daydiff60 = [];
nightdiff60 = [];
bbpavg60 = repmat(bbpttl60, 1, 12);
for i = 1:length(dayplot60)
    daydiff60 = [daydiff60 (dayplot60(i) - bbpavg60(i))];
    nightdiff60 = [nightdiff60 (nightplot60(i) - bbpavg60(i))];
end

fig8 = figure(8);
hold on
plot([1 166], [0 0],'color',[.85 .85 .85],'linewidth',2,'HandleVisibility','off');
plot(71:166, floatdiff60,'c','linewidth',1);
plot(daydiff60,'--b','linewidth',2);
plot(nightdiff60,':r','linewidth',2);
axis tight
hold off
xticks([1:12:166]);
xlab = 2006:2019;
xticklabels(xlab);
legend('Float','Satellite (day)','Satellite (night)');
xlabel('date');
ylabel('deviation from monthly mean (m^{-1})');
title('bbp anomaly 60\circ (2012-2019)');
print(fig8, '/Users/aamatya/Desktop/Summer 2019/8-anomaly60','-dpdf');

%% Monthly satellite bbp anomalies - by latitude 
%6/2006 to 4/2017 

fig9 = figure(9);
hold on
line([1 137],[0 0],'Color',[.9 .9 .9],'linewidth',5,'Handlevisibility','off');
plot(daydiff,'Color',[.8 .4 .4] ,'linewidth',1);
plot(daydiff60,'Color',[.4 .4 .8],'linewidth',1); %length w/o nans = 131
legend('b_{bp} Anomaly, below 25\circ S', 'b_{bp} Anomaly, below 60\circ S');
title('Monthly LIDAR b_{bp} Anomaly, 25\circS and 60\circS (2006-2017)')
xlim([5 137]);
ylim([-.0005 .00085]);
xticks([1:12:166]); xlab = 2006:2019; xticklabels(xlab);
xlabel('Year');
ylabel('Deviation from Monthly b_{bp} Average, 2006-2019  (m^{-1})');
hold off
print(fig9, '/Users/aamatya/Desktop/Summer 2019/9-anomalybothlats','-dpdf');


