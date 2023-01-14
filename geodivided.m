%% Climatology based on region
% Pacific: [125 to 180], [-180 to -70]
% Indian: [25 to 125]
% Atlantic: [-70 to 25]
%% Load Data

cd('/Users/aamatya/Desktop/SOCCOM');
load Datosint2.mat;     % Float bbp
load Datoslid_SO.mat;   % Satellite bbp
load sst_anomalies.mat; % Temperature

%Simplify time variables into single-rowed arrays
for i = 1:length(Datosint2)
    Datosint2(i).month_simple = Datosint2(i).month(1,:);
    Datosint2(i).year_simple = Datosint2(i).year(1,:);
end

%Extract desired variables
satmonth = extractfield(Datoslid_SO,'month_simple');
satbbp = extractfield(Datoslid_SO, 'bbp');
sattime = extractfield(Datoslid_SO, 'time_combined');
satmonthrep = month(extractfield(Datoslid_SO, 'time'));
satyearrep = year(extractfield(Datoslid_SO, 'time'));
satlat = extractfield(Datoslid_SO, 'lat');
satlon = extractfield(Datoslid_SO, 'lon');
satphyto = extractfield(Datoslid_SO, 'cphyto');
%% indices
% 1 - Between 30 and 55, Pacific
geodiv(1).sidx = find((satlat >= -55) & (satlat <= -30) & ((satlon <= -70) | (satlon >= 125)));
% 2 - Between 30 and 55, Indian
geodiv(2).sidx = find((satlat >= -55) & (satlat <= -30) & (satlon >= 25) & (satlon <= 125));
% 3 - Between 30 and 55, Atlantic
geodiv(3).sidx = find((satlat >= -55) & (satlat <= -30) & (satlon >= -70) & (satlon <= 25));
% 4 - Below 55, Pacific
geodiv(4).sidx = find((satlat <= -55) & ((satlon <= -70) | (satlon >= 125)));
% 5 - Below 55, Indian
geodiv(5).sidx = find((satlat <= -55) & (satlon >= 25) & (satlon <= 125));
% 6 - Below 55, Atlantic
geodiv(6).sidx = find((satlat <= -55) & (satlon >= -70) & (satlon <= 25));
% 7 - All, Pacific
geodiv(7).sidx = find((satlon <= -70) | (satlon >= 125));
% 8 - All, Indian
geodiv(8).sidx = find((satlon >= 25) & (satlon <= 125));
% 9 - All, Atlantic
geodiv(9).sidx = find((satlon >= -70) & (satlon <= -25));
%10 - Between 30 and 55, All
geodiv(10).sidx = find((satlat >= -55) & (satlat <= -30));
%11 - Below 55, All
geodiv(11).sidx = find(satlat <= -55);
%12 - All, All
geodiv(12).sidx = 1:length(satlat);

for j = 1:12 % Each index
    geodiv(j).month = satmonthrep(geodiv(j).sidx);
    geodiv(j).year = satyearrep(geodiv(j).sidx);
    geodiv(j).bbp = satbbp(geodiv(j).sidx);
    %geodiv(j).time = sattime(geodiv(j).sidx);
    geodiv(j).cphyto = satphyto(geodiv(j).sidx);
end

for i = 1:12 % Each index
    for j = 1:12 % Each month
        b = find(geodiv(i).month == j);
        geodiv(i).bbpmonthly(j).month = geodiv(i).bbp(b);
        geodiv(i).phytomonthly(j).month = geodiv(i).cphyto(b);
        clear b;
        geodiv(i).bbpmonthly(j).average = nanmean(geodiv(i).bbpmonthly(j).month);
        geodiv(i).phytomonthly(j).average = nanmean(geodiv(i).phytomonthly(j).month);
    end
end

for i = 1:12 % Each index
    tempmonth = geodiv(i).month;
    tempyear = geodiv(i).year;
    tempbbp = geodiv(i).bbp;
    tempphyto = geodiv(i).cphyto;
    for j = 1:12 %month
        for k = 1:12 %year
            thisdate = find((tempmonth == j) & (tempyear - 2005 == k));
            geodiv(i).bbpmonthly(j).peryear(k) = nanmean(tempbbp(thisdate));
            geodiv(i).phytomonthly(j).peryear(k) = nanmean(tempphyto(thisdate));
        end
    end
end

for i = 1:12 % Each index
    geodiv(i).plot = [];
    geodiv(i).phytoplot = [];
    for j = 1:12 %month
        for k = 1:12 %year
            thisday = geodiv(i).bbpmonthly(k).peryear(j);
            thisdayphyto = geodiv(i).phytomonthly(k).peryear(j);
            geodiv(i).plot = [geodiv(i).plot thisday];
            geodiv(i).phytoplot = [geodiv(i).phytoplot thisdayphyto];
        end
    end
end
clear thisdate thisday thisnight day night tempmonth tempyear tempbbp thisdayphyto;

for i = 1:12 % Each index
    avgarray = [];
    avgarrayphyto = [];
    for j = 1:12
        avgarray = [avgarray geodiv(i).bbpmonthly(j).average];
        avgarrayphyto = [avgarrayphyto geodiv(i).phytomonthly(j).average];
    end
    geodiv(i).avgarray = avgarray
    geodiv(i).avgarrayphyto = avgarrayphyto
end

for i = 1:12 % Each index
    geodiff = [];
    geodiffphyto = [];
    monthavg = repmat(geodiv(i).avgarray,1,12);
    monthavgphyto = repmat(geodiv(i).avgarrayphyto, 1, 12);
    for k = 1:length(geodiv(i).plot)
        geodiff = [geodiff (geodiv(i).plot(k) - monthavg(k))];
        geodiffphyto = [geodiffphyto (geodiv(i).phytoplot(k) - monthavgphyto(k))];
    end
    geodiv(i).geodiff = geodiff; 
    geodiv(i).geodiffphyto = geodiffphyto;
end

%% Average Plot
%addpath('/Users/aamatya/Desktop/suplabel');
fig12 = figure(12);
thetitles = {'30-55\circ, Pacific','30-55\circ, Indian','30-55\circ, Atlantic','55\circ, Pacific'...
    '55\circ, Indian','55\circ, Atlantic','All, Pacific','All, Indian','All, Atlantic'};

for i = 1:9 % Plot the first 9 indices
    ax(i) = subplot(3,3,i);
    plot(geodiv(i).plot);
    xlim([5 137]);
    xticks([1:24:144+12]); xlab = (2006:2:2018); xticklabels(xlab);
    ylim([.0013 .0035]);
end

titles = {'Pacific','Indian','Atlantic','30-55\circ','55\circ','All'};
for i = 1:3;
    subplot(ax(i));
    title(titles(i),'fontsize',14,'fontweight','normal','color',[.8 .3 .3],'fontname','times');
end

subplot(ax(1));
y = ylabel('30\circS-55\circS','fontsize',14,'color',[.3 .7 .4],'fontname','times'); 
set(get(gca,'ylabel'),'rotation',0); set(y, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
subplot(ax(4));
y = ylabel('<55\circS','fontsize',14,'color',[.3 .7 .4],'fontname','times');
set(get(gca,'ylabel'),'rotation',0); set(y, 'Units', 'Normalized', 'Position', [-0.16, 0.5, 0]);
subplot(ax(7));
y = ylabel('All','fontsize',14,'color',[.3 .7 .4],'fontname','times');
set(get(gca,'ylabel'),'rotation',0); set(y, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
    
[ax, h1] = suplabel('Monthly LIDAR b_{bp} Averages by Lat and Lon (2006 - 2017)','t'); 
set(h1,'fontsize',17,'fontweight','normal','fontname','times');
[ax, h2] = suplabel('Year','x'); 
set(h2,'fontsize',17,'fontweight','normal','fontname','times');
[ax, h3] = suplabel('Monthly b_{bp} Average (m^{-1})','y'); 
set(h3,'fontsize',17,'fontweight','normal','fontname','times');

set(fig12,'PaperUnits', 'centimeters', 'PaperSize', [50 25], 'PaperPosition', [0 0 50 25])
print(fig12, '/Users/aamatya/Desktop/Summer 2019/12-geodiv','-dpdf');
%% Anomaly Plot
%color back by qudrant
fig13 = figure(13);

for i = 1:9 % Plot the first nine indices 
    ax(i) = subplot(3,3,i)
    hold on
    plot([1 144],[0 0],'color',[.93 .93 .93],'linewidth',4);
    plot(geodiv(i).geodiff);
    xticks([1:24:144+12]); xlab = (2006:2:2018); xticklabels(xlab);
    ylim([-.0006 .0008]);
    hold off
end

titles = {'Pacific','Indian','Atlantic','30-55\circ','55\circ','All'};
for i = 1:3;
    subplot(ax(i));
    title(titles(i),'fontsize',14,'fontweight','normal','color',[.8 .3 .3],'fontname','times');
end

subplot(ax(1));
y = ylabel('30\circS-55\circS','fontsize',14,'color',[.3 .7 .4],'fontname','times'); 
set(get(gca,'ylabel'),'rotation',0); set(y, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
subplot(ax(4));
y = ylabel('<55\circS','fontsize',14,'color',[.3 .7 .4],'fontname','times');
set(get(gca,'ylabel'),'rotation',0); set(y, 'Units', 'Normalized', 'Position', [-0.16, 0.5, 0]);
subplot(ax(7));
y = ylabel('All','fontsize',14,'color',[.3 .7 .4],'fontname','times');
set(get(gca,'ylabel'),'rotation',0); set(y, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
    
[ax, h1] = suplabel('Monthly LIDAR b_{bp} Anomalies by Lat and Lon (2006 - 2017)','t'); 
set(h1,'fontsize',17,'fontweight','normal','fontname','times');
[ax, h2] = suplabel('Year','x'); 
set(h2,'fontsize',17,'fontweight','normal','fontname','times');
[ax, h3] = suplabel('b_{bp} Anomaly (m^{-1})','y'); 
set(h3,'fontsize',17,'fontweight','normal','fontname','times');

set(fig13,'PaperUnits', 'centimeters', 'PaperSize', [50 25], 'PaperPosition', [0 0 50 25])
print(fig13, '/Users/aamatya/Desktop/Summer 2019/13-geoanomaly','-dpdf');



