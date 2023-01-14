%% BBP Anomaly
% Simplify time variables into single-rowed arrays
for i = 1:length(Datosint2)
    Datosint2(i).month_simple = Datosint2(i).month(1,:);
    Datosint2(i).year_simple = Datosint2(i).year(1,:);
end

% Extract desired variables
floatyear = extractfield(Datosint2, 'year_simple');
b = find(floatyear >= 2012 & floatyear <= 2017);
floatmonth = extractfield(Datosint2,'month_simple');
floatmonth = floatmonth(b);
floatbbp = extractfield(Datosint2,'bbp_mld');
floatbbp = floatbbp(b);
floatlon = extractfield(Datosint2,'lon');
floatlon = floatlon(b);
floatqf = extractfield(Datosint2, 'lat_QF');
floatqf = floatqf(b);
floatlat = extractfield(Datosint2, 'lat');
floatlat = floatlat(b);

satyearrep = year(extractfield(Datoslid_SO, 'time'));
a = find(satyearrep >= 2012 & satyearrep <= 2017);
satbbp = extractfield(Datoslid_SO, 'bbp');
satbbp = satbbp(a);
satmonthrep = month(extractfield(Datoslid_SO, 'time'));
satmonthrep = satmonthrep(a);
satlat = extractfield(Datoslid_SO, 'lat');
satlat = satlat(a);
satlon = extractfield(Datoslid_SO, 'lon');
satlon = satlon(a);
b = find(satlon <= 0);
satlon(b) = satlon(b) + 360;
%% 2012 - 2017

for i = 1:12
    a = find(satmonthrep == i);
    bbp(i).clim(1:360, 1:65) = NaN;
    bbp(i).total = satbbp(a);
    bbp(i).lat = satlat(a);
    bbp(i).lon = satlon(a);
end
%%
lonprev = 1; %ceil(bbp(1).lon(1));
latprev = 1; %floor(bbp(1).lat(1))+90;
bbp(1).clim(lonprev, latprev) = NaN; %bbp(1).total(1);
count = 1;
for i = 1:12
    disp(i)
    for j = 1:length(bbp(i).lon)
        thelat = round(bbp(i).lat(j)) + 90;
        thelon = round(bbp(i).lon(j));
        if thelon >= 360
            thelon = 360;
        elseif thelon == 0
            thelon = 1;
        end
        
        if (thelon == lonprev & thelat == latprev) 
            count = count+1;
        else 
            bbp(i).clim(lonprev, latprev) = bbp(i).clim(lonprev, latprev) / count;
            bbp(i).clim(thelon, thelat) = 0;
            count = 1;
        end
        
        bbp(i).clim(thelon, thelat) = bbp(i).clim(thelon, thelat) + (bbp(i).total(j));
        lonprev = thelon;
        latprev = thelat;
    end
end
%%

maplats_bbp = [];
maplons_bbp = [];
mapanoms_bbp = [];

for i = 1:length(Datosint2)
    for j = 1:length(Datosint2(i).lat)
        
        thelat = round(Datosint2(i).lat(j)) + 90;
        thelon = round(Datosint2(i).lon(j));
        thetime = Datosint2(i).month_simple(j);
        
        if isnan(thelat) == 1 || isnan(thelon) == 1
            continue
        elseif thelon == 0
            thelon = 1;
        elseif thelon >= 360
            thelon = 360;
        end
        
        thebase = bbp(thetime).clim(thelon, thelat);
        
        if (isnan(thebase) == 1)
            continue
        end
        
        maplats_bbp = [maplats_bbp Datosint2(i).lat(j)];
        maplons_bbp = [maplons_bbp Datosint2(i).lon(j)];
        mapanoms_bbp = [mapanoms_bbp Datosint2(i).bbp_mld(j) - thebase];
    end
end

filename = gunzip('gshhs_c.b.gz', tempdir);
shorelines = gshhs(filename{1});
delete(filename{1})
levels = [shorelines.Level];
land = (levels == 1);

fig18 = figure(18);
worldmap([-90 -25],[0 360])
scatterm(maplats_bbp , maplons_bbp,[], mapanoms_bbp,'filled'); colormap(redblue),colorbar
geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
title('bbp Anomalies by float')
caxis([-0.0015 0.0015]);

% print(fig18, '/Users/aamatya/Desktop/Summer 2019/18-bbpanomalymap', '-dpdf');
%% Float by year
maplats_bbp = [];
maplons_bbp = [];
mapanoms_bbp = [];

for i = 1:6 %years
    for j = 1:length(Datosint2)
        for k = 1:length(Datosint2(j).lat)
            if year(Datosint2(j).month_simple(k)) - 2011 ~= i
                continue
            end
            
            thelat = floor(Datosint2(j).lat(k)) + 90;
            thelon = ceil(Datosint2(j).lon(k));
            thetime = Datosint2(j).month_simple(k);
            
            if isnan(thelat) | isnan(thelon)
                continue
            elseif thelon == 0
                thelon = 1;
            elseif thelon >= 360
                thelon = 360;
            end
            
            thebase = bbp(thetime).clim(thelon, thelat);
            
            if isnan(thebase);
                continue
            end
            
            maplats_bbp = [maplats_bbp (Datosint2(j).lat(k))];
            maplons_bbp = [maplons_bbp Datosint2(j).lon(k)];
            mapanoms_bbp = [mapanoms_bbp Datosint2(j).bbp(k) - thebase];
        end
    end
    yearlybbp(i).lat = maplats_bbp;
    yearlybbp(i).lon = maplons_bbp;
    yearlybbp(i).anom = mapanoms_bbp;
end


fig20 = figure(20);
for i = 1:6
    subplot(2,3,i);
    worldmap([-90 -25],[0 360])
    scatterm(yearlybbp(i).lat , yearlybbp(i).lon,[], yearlybbp(i).anom,'filled'); colormap(redblue),colorbar
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
    title('bbp Anomalies by float')
    caxis([-0.002 0.002]);
end

%% Satellite all
for k = 1:12
    maplats_bbpsat = [];
    maplons_bbpsat = [];
    mapanoms_bbpsat = [];
    maplats_bbpsat_tot = [];
    maplons_bbpsat_tot = [];
    mapanoms_bbpsat_tot = [];
    disp(k)
    for i = 1:length(Datoslid_SO)
        id = round(rand(1000, 1) * length(Datoslid_SO(i).lat));
        for j = 1:1000 %:length(Datoslid_SO(i).lat)
            if id(j) == 0
                continue
            elseif month(Datoslid_SO(i).time(id(j))) ~= k
                continue
            end
            thelat = floor(Datoslid_SO(i).lat(id(j))) + 90;
            thelon = ceil(Datoslid_SO(i).lon(id(j)));
            thetime = month(Datoslid_SO(i).time(id(j)));
            
            if isnan(thelat) | isnan(thelon)
                continue
            elseif thelon <= 0
                thelon = thelon + 360;
            end
            
            thebase = bbp(thetime).clim(thelon, thelat);
            
            if isnan(thebase);
                continue
            end
            
            maplats_bbpsat = [maplats_bbpsat Datoslid_SO(i).lat(id(j))];
            maplons_bbpsat = [maplons_bbpsat Datoslid_SO(i).lon(id(j))];
            mapanoms_bbpsat = [mapanoms_bbpsat Datoslid_SO(i).bbp(id(j)) - thebase];
            maplats_bbpsat_tot = [maplats_bbpsat_tot Datoslid_SO(i).lat(id(j))];
            maplons_bbpsat_tot = [maplons_bbpsat_tot Datoslid_SO(i).lon(id(j))];
            mapanoms_bbpsat_tot = [mapanoms_bbpsat_tot Datoslid_SO(i).bbp(id(j))];
        end
    end
    yearlybbp(k).satlat = maplats_bbpsat;
    yearlybbp(k).satlon = maplons_bbpsat;
    yearlybbp(k).satanom = mapanoms_bbpsat; 
    yearlybbp(k).satlat_tot = maplats_bbpsat_tot;
    yearlybbp(k).satlon_tot = maplons_bbpsat_tot;
    yearlybbp(k).satanom_tot = mapanoms_bbpsat_tot;
end
% 
fig21 = figure(21);
for k = 1:12
    subplot(3,4,k)
    worldmap([-90 -25],[0 360])
    scatterm(yearlybbp(k).satlat, yearlybbp(k).satlon, [], yearlybbp(k).satanom,'filled'); colormap(redblue),colorbar
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
    title('bbp satellite anomalies')
    caxis([-0.0006 0.0006]);
end


fig22 = figure(22);
for k = 1:12
    subplot(3,4,k)
worldmap([-90 -25],[0 360])
scatterm(yearlybbp(k).satlat_tot, yearlybbp(k).satlon_tot, [], yearlybbp(k).satanom_tot,'filled'); colormap(summer),colorbar
geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
title('bbp total monthly satellite')
caxis([0 .004]); 
end
print(fig22, '/Users/aamatya/Desktop/Summer 2019/22 - bbp total sat', '-dpdf');

% figure,
% % for i = 1:6
% %     subplot(2,3,i)
% worldmap([-90 -25],[0 360])
% scatterm(maplats_bbpsat_tot, maplons_bbpsat_tot, [], mapanoms_bbpsat_tot,'filled'); colormap(bone),colorbar
% geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
% title('bbp Anomalies by float')
% caxis([0.0025 .004]);
% % end
%% Float nonanomaly
figure,
for k = 1:6
    maplats_bbpsat_tot = [];
    maplons_bbpsat_tot = [];
    mapanoms_bbpsat_tot = [];
    for i = 1:length(Datoslid_SO)
        disp(i)
        id = round(rand(1000, 1) * length(Datoslid_SO(i).lat));
        for j = 1:1000 %:length(Datoslid_SO(i).lat)
            if year(Datoslid_SO(i).time(id(j))) - 2011 ~= k
                continue
            end
            
            maplats_bbpsat_tot = [maplats_bbpsat_tot Datoslid_SO(i).lat(id(j))];
            maplons_bbpsat_tot = [maplons_bbpsat_tot Datoslid_SO(i).lon(id(j))];
            mapanoms_bbpsat_tot = [mapanoms_bbpsat_tot Datoslid_SO(i).bbp(id(j))];
        end
    end
    yearlybbp(k).satlat_tot = maplats_bbpsat_tot;
    yearlybbp(k).satlon_tot = maplons_bbpsat_tot;
    yearlybbp(k).satanom_tot = mapanoms_bbpsat_tot;
end

for i = 1:6
    subplot(2,3,i)
    worldmap([-90 -25],[0 360])
    scatterm(maplats_bbpsat_tot, maplons_bbpsat_tot, [], mapanoms_bbpsat_tot,'filled'); colormap(redblue),colorbar
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
    title('bbp Anomalies by float')
    caxis([0 .005]);
end