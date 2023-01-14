%% Nitrate Anomaly (March 2012 - March 2019)

cd('/Users/aamatya/Desktop/SOCCOM/nitdata');
list = dir('/Users/aamatya/Desktop/SOCCOM/nitdata');

nitdata(1).idx = 35:60; % between -30 and -55
nitdata(2).idx = 1:35; % below -55
nitdata(3).idx = 1:60; % all below -30

for i = 1:3
    for j = 4:length(list) %each month
        nitrate_mean = ncread(list(j).name,'n_an');
        nitrate_mean = nitrate_mean(:,nitdata(i).idx,1:11); %below 30S
        fillval_mean = ncreadatt(list(j).name,'n_an','_FillValue');
        nitrate_mean(nitrate_mean == fillval_mean) = NaN;
        nitdata(i).overall(j-3) = nanmean(nitrate_mean(:));
    end
end
nitdata = rmfield(nitdata,'idx'); clear i j fillval_mean list nitrate_mean;

for i = 1:length(Datosint2)
    Datosint2(i).month_simple = Datosint2(i).month(1,:);
    Datosint2(i).year_simple = Datosint2(i).year(1,:);
end
floatmonth = extractfield(Datosint2,'month_simple');
floatyear = extractfield(Datosint2, 'year_simple');
floatnitrate = extractfield(Datosint2, 'nitrate_mld');
floatlatt = extractfield(Datosint2, 'lat');
floatqf = extractfield(Datosint2, 'lat_QF');

floatlat = floatlatt(floatqf  ~= 1);
nitdata(1).fltidx = find((floatlat <= -30) & (floatlat >=-55));
nitdata(2).fltidx = find(floatlat <= -55);
nitdata(3).fltidx = 1:length(floatlat);

for i = 1:3
    nitdata(i).floatmonth = floatmonth(nitdata(i).fltidx);
    nitdata(i).floatyear = floatyear(nitdata(i).fltidx);
    nitdata(i).floatnitrate = floatnitrate(nitdata(i).fltidx);
end

for i = 1:3
    for j = 1:12
        b = find(nitdata(i).floatmonth == j);
        nitdata(i).nitmonthly(j).month = nitdata(i).floatnitrate(b);
        clear b
        nitdata(i).nitmonthly(j).average = nanmean(nitdata(i).nitmonthly(j).month);
    end
end

for i = 1:3
    tempmonth = nitdata(i).floatmonth;
    tempyear = nitdata(i).floatyear;
    tempnitrate = nitdata(i).floatnitrate;
    for j = 1:12 %month
        display(j)
        for k = 1:8 %year
            thisdate = find((tempmonth == j) & (tempyear - 2011 == k));
            nitdata(i).nitmonthly(j).peryear(k) = nanmean(tempnitrate(thisdate));
        end
    end
end

for i = 1:3
    nitdata(i).plot = [];
    for j = 1:8 %year
        for k = 1:12 %month
            thisday = nitdata(i).nitmonthly(k).peryear(j);
            nitdata(i).plot = [nitdata(i).plot thisday];
        end
    end
end
clear thisdate thisday thisnight day night tempmonth tempyear tempbbp tempnitrate ;

for i = 1:3
    nit_base = repmat(nitdata(i).overall,1,8);
    nitdiff = [];
    for k = 1:length(nitdata(i).plot)
        nitdiff = [nitdiff (nitdata(i).plot(k) - nit_base(k))];
    end
    nitdata(i).nitdiff = nitdiff;
end

nitdata = rmfield(nitdata, 'fltidx');
%% Anomaly plot (nitrate, temperature, bbp) per latitude
% Latitudes = -55 to -30, below -55, all 

% load sst_anomalies.mat;
firstidx_f = find(dates == datenum('01-Jan-2012')); %jun 2006
sstdata_f = data(firstidx_f:length(dates)); %until 10/18

fig15 = figure(15); %start 3/12
for i = 1:3 % Per latitude
    sb(i) = subplot(1,3,i);
    hold on
    nitdiff = nitdata(i).nitdiff;
    
    yyaxis right
%   Fix legend to add temp if following line is uncommented
%   plot(1:length(sstdata_f), sstdata_f,'-','color',[.3 .7 .4]); 
    plot(1:70, geodiv(i+9).geodiff(75:144).*2000);
    ylim([-1 1.5]);
    
    yyaxis left
    plot(1:length(nitdiff),nitdiff,'-b');
    ylim([-10 15]);
    
    xlim([5 96]);
    xticks([1:12:96]); xlab = 2012:2019; xticklabels(xlab);
    xlabel('Year');
    
    legend('bbp','nitrate');
    titles = {'between -55 -30','below -55','all latitudes'};
    title(titles(i));
    hold off
    
end


set(fig15,'PaperUnits', 'centimeters', 'PaperSize', [50 25], 'PaperPosition', [0 0 50 25])
print(fig15, '/Users/aamatya/Desktop/Summer 2019/15-nitrateanomaly','-dpdf');

clear firstidx_f i j k sb;

%% Nitrate Anomaly Map

for j = 4:length(list) % each month
    nitrate_mean = ncread(list(j).name,'n_an');
    nitrate_mean = nitrate_mean(:,1:65,1:11); %below 30S
    fillval_mean = ncreadatt(list(j).name,'n_an','_FillValue');
    nitrate_mean(nitrate_mean == fillval_mean) = NaN;
    nitdata(4).overall(j-3).grid = nanmean(nitrate_mean,3); % Below 25S
end

for i = 1:length(Datosint2)
    Datosint2(i).time_combined = Datosint2(i).time(1,:);
    Datosint2(i).lat(Datosint2(i).lat_QF == 1) = NaN;
    Datosint2(i).lon(Datosint2(i).lat_QF == 1) = NaN;
    Datosint2(i).nitrate_mld(Datosint2(i).lat_QF == 1) = NaN;
end

maplats = [];
maplons = [];
mapanoms = [];
maptimes = [];

for i = 1:length(Datosint2)
    for j = 1:length(Datosint2(i).lat)
        thelat = round(Datosint2(i).lat(j)) + 90;
        thelon = round(Datosint2(i).lon(j));
        thetime = month(Datosint2(i).time_combined(j));
        
        if isnan(thelat) | isnan(thelon)
            continue
        elseif thelon == 0
            thelon = 1;
        elseif thelon >= 360
            thelon = 360;
        end
        
        thebase = nitdata(4).overall(thetime).grid(thelon, thelat);
        
        if isnan(thebase);
            continue
        end
        
        maplats = [maplats Datosint2(i).lat(j)];
        maplons = [maplons Datosint2(i).lon(j)];
        mapanoms = [mapanoms Datosint2(i).nitrate_mld(j) - thebase];
        maptimes = [maptimes Datosint2(i).time_combined(j)];
    end
end

filename = gunzip('gshhs_c.b.gz', tempdir);
shorelines = gshhs(filename{1});
delete(filename{1})
levels = [shorelines.Level];
land = (levels == 1);

fig16 = figure(16)
worldmap([-90 -25],[0 360])
scatterm(maplats,maplons,[], mapanoms,'filled'); colormap(redblue),colorbar
geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
title('Nitrate Anomalies by float')
caxis([-10 10]);

% print(fig16, '/Users/aamatya/Desktop/Summer 2019/16-nitrateanomalymap','-dpdf');
%% Float nitrate anomaly per month
for k = 1:12
    maplats_mon = [];
    maplons_mon = [];
    mapanoms_mon = [];
    for i = 1:length(Datosint2)
        for j = 1:length(Datosint2(i).lat)
            if Datosint2(i).month_simple(j) ~= k
                continue
            end
            thelat = round(Datosint2(i).lat(j)) + 90;
            thelon = round(Datosint2(i).lon(j));
            thetime = month(Datosint2(i).time_combined(j));
            
            if isnan(thelat) | isnan(thelon)
                continue
            elseif thelon == 0
                thelon = 1;
            elseif thelon >= 360
                thelon = 360;
            end
            
            thebase = nitdata(4).overall(thetime).grid(thelon, thelat);
            
            if isnan(thebase);
                continue
            end
            
            maplats_mon = [maplats_mon Datosint2(i).lat(j)];
            maplons_mon = [maplons_mon Datosint2(i).lon(j)];
            mapanoms_mon = [mapanoms_mon Datosint2(i).nitrate_mld(j) - thebase];
        end
    end
    nitdata(4).monthlymap(k).lat = maplats_mon;
    nitdata(4).monthlymap(k).lon = maplons_mon;
    nitdata(4).monthlymap(k).anom = mapanoms_mon;
end

% fig23 =figure(23);
% nc
% print(fig23, '/Users/aamatya/Desktop/Summer 2019/23 - nitrate anom monthly', '-dpdf');
%% climatology per month

for i = 1:12
    for j = 1:360
        for k = 1:65
            nitdata(i).total = [nitdata(i).total nitdata(4).overall(i).grid(j,k)];
        end
    end
end

trashlat = [1:65 1:65]
for i = 1:12
    subplot(3,4,i);
    
    worldmap([-90 -25],[0 360])
    scatterm(repmat(1:65,1, 360),repmat(1:360, 1, 65),[], nitdata(4).overall(i).grid,'filled'); colormap(redblue),colorbar
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
    title('Nitrate Anomalies by float')
    caxis([-10 10]);
end








