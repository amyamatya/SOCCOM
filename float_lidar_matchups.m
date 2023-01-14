function [Datosint2] = float_lidar_matchups(Datoslid_SO,Datosint)
cd('/Users/aamatya/Desktop/SOCCOM');
load Datos_interpolated_12Mar2019.mat
load Datos_lidar.mat
for n = 1:length(Datoslid_SO)
    if ~isempty(Datoslid_SO(n).time)
        Datoslid_SO(n).year_simple = year(Datoslid_SO(n).time(1));
        Datoslid_SO(n).month_simple = month(Datoslid_SO(n).time(1));
        if ~isempty(Datoslid_SO(n).year_simple)
            Datoslid_SO(n).time_combined = datenum(datetime(Datoslid_SO(n).year_simple,Datoslid_SO(n).month_simple,01));
        else
            Datoslid_SO(n).time_combined = NaN;
        end
    else
        Datoslid_SO(n).year_simple = NaN;
        Datoslid_SO(n).month_simple = NaN;
    end
end
clear n;

bbp700 = extractfield(Datoslid_SO, 'bbp532').* (532./700);

for n = 1:length(Datoslid_SO)
   if isempty(Datoslid_SO(n).year_simple) == 0
       Datoslid_SO(n).time_combined = datenum(datetime(Datoslid_SO(n).year_simple,Datoslid_SO(n).month_simple,01));
   else
       Datoslid_SO(n).time_combined = NaN;
   end
end

%% test
sat_times = extractfield(Datoslid_SO,'time_combined');
for ii = 1:length(Datosint)
    display(ii)
    [a,b] = size(Datosint(ii).time);
    for n = 1:b
        %display(n)
        float_month = month(Datosint(ii).time(1,n));
        float_year = year(Datosint(ii).time(1,n));
        float_day = day(Datosint(ii).time(1,n));
        float_lat = Datosint(ii).lat(1,n);
        float_lon = Datosint(ii).lon(1,n);
        float_latQF = Datosint(ii).lat_QF(1,n);
        if float_latQF == 1
            float_lat = NaN;
            float_lon = NaN;
        end
        if float_lon > 180
            float_lon = float_lon - 360;
        end
        idx = find(month(sat_times) == float_month & year(sat_times) == float_year);
        if isempty(idx) == 0
            idx = idx(1);
            sat_date = Datoslid_SO(idx).time;
            sat_day = day(Datoslid_SO(idx).time);
            day_diff = float_day - sat_day;
            same_day_idx = find(day_diff == 0);
            if isempty(same_day_idx) == 0
                thedays = sat_date(same_day_idx);
                sat_lats = Datoslid_SO(idx).lat(same_day_idx);
                sat_lons = Datoslid_SO(idx).lon(same_day_idx);
                sat_bbp700 = Datoslid_SO(idx).bbp532(same_day_idx) .* (532./700);
                latlon1 = [float_lat float_lon];
                for jj = 1:length(thedays)
                    latlon2 = [sat_lats(jj) sat_lons(jj)];
                    [d1km d2km] = lldistkm(latlon1,latlon2);
                    thedistance(jj) = d1km;
                end
                outof_distance = find(thedistance > 200);
                %isempty(outof_distance)
                meand(n) = nanmean(thedistance);
                sat_lats(outof_distance) = NaN;
                sat_lons(outof_distance) = NaN;
                sat_bbp700(outof_distance) = NaN;
                thedays(outof_distance) = NaN;
                sat_lats_final(n) = nanmean(sat_lats);
                sat_lons_final(n) = nanmean(sat_lons);
                sat_bbp700_final(n) = nanmean(sat_bbp700);
                thedays_final(n) = nanmean(thedays);
            else
                sat_lats_final(n) = NaN;
                sat_lons_final(n) = NaN;
                sat_bbp700_final(n) = NaN;
                thedays_final(n) = NaN;
                meand(n) = NaN;
            end
        else
            sat_lats_final(n) = NaN;
            sat_lons_final(n) = NaN;
            sat_bbp700_final(n) = NaN;
            thedays_final(n) = NaN;
            meand(n) = NaN;
        end
        clearvars -except sat_lats_final sat_lons_final sat_bbp700_final thedays_final meand ii Datoslid_SO Datosint sat_times
    end
    %newidx = find(nonzeros(sat_bbp700_final));
    Datosint(ii).lat_sat_200 = sat_lats_final;%(newidx);
    Datosint(ii).lon_sat_200 = sat_lons_final;%(newidx);
    Datosint(ii).bbp700_sat_200 = sat_bbp700_final;%(newidx);
    Datosint(ii).date_sat_200 = thedays_final;%(newidx);
    Datosint(ii).meand_sat_200 = meand;%(newidx);
    clearvars -except Datoslid_SO Datosint sat_times 
end
Datosint2 = Datosint;
return
for jk = 1:length(Datosint)
    [sizt, sizy] = size(Datosint(jk).time);
    for kk = 1:sizy
        themean = nanmean(Datosint(jk).meand_sat_200(1:kk));
        Datosint(jk).meandist = themean;
    end
end

return

%% Analysis [Figures]
meandist1 = extractfield(Datosint2, 'meand_sat');
lat1 = extractfield(Datosint2, 'lat_sat');
bbpsat1 = extractfield(Datosint2,'bbp700_sat');
bbp1 = extractfield(Datosint2,'bbp_mld');
lon1 = extractfield(Datosint2, 'lon_sat');
bbp1 = extractfield(Datosint2,'bbp_mld');
date1 = extractfield(Datosint2,'date_sat');
month1 = datetime(date1, 'ConvertFrom','datenum');
month1 = month(month1);
%%
fig1 = figure(1);
scatter(log(bbp1),log(bbpsat1))
title('satellite vs. float bbp');
set(fig1,'PaperUnits', 'centimeters', 'PaperSize', [20 10], 'PaperPosition', [0 0 20 10])
print(fig1, 'Figure 1', '-dpdf');

fig2 = figure(2);
summer1 = find((month1 >= 9 & month1 <= 12) | (month1 >= 1 & month1 <= 3));
winter1 = find(month1 > 3 & month1 < 9);
winterbbp1 = bbpsat1(winter1);
summerbbp1 = bbpsat1(summer1);
winterbbp2 = bbp1(winter1);
summerbbp2 = bbp1(summer1);
hold on
scatter(log(winterbbp2), log(winterbbp1),'r');
scatter(log(summerbbp2), log(summerbbp1),'b');
legend('Winter', 'Summer');
title('by season');
hold off
set(fig2,'PaperUnits', 'centimeters', 'PaperSize', [20 10], 'PaperPosition', [0 0 20 10])
print(fig2, 'Figure 2', '-dpdf');

fig3 = figure(3);
hold on
view(3);
scatter3(lat1, lon1, bbp1);
scatter3(lat1, lon1, bbpsat1);
legend('float','satellite');
rotate3d on
xlabel('latitude');
ylabel('longitude');
zlabel('bbp');
grid on
hold off
set(fig3,'PaperUnits', 'centimeters', 'PaperSize', [20 10], 'PaperPosition', [0 0 20 10])
print(fig3, 'Figure 3', '-dpdf');

for ii = 1:length(Datosint)
Datosint(ii).time_simple = Datosint(ii).time(1,:);
end
date2 = extractfield(Datosint,'time_simple');

fig4 = figure(4);
scatter(log(bbp1),log(bbpsat1),[],month(date2))
colorbar
colormap(jet)
set(fig4,'PaperUnits', 'centimeters', 'PaperSize', [20 10], 'PaperPosition', [0 0 20 10])
print(fig4, 'Figure 4', '-dpdf');

fig5 = figure(5);
scatter(log(bbp1),log(bbpsat1),[],lat1);
title('by latitude');
colorbar;
colormap(jet);
set(fig5,'PaperUnits', 'centimeters', 'PaperSize', [20 10], 'PaperPosition', [0 0 20 10])
print(fig5, 'Figure 5', '-dpdf');

fig6 = figure(6);
scatter(log(bbp1),log(bbpsat1),[],lon1);
title('by longitude');
colorbar;
colormap(jet);
%set(fig6,'PaperUnits', 'centimeters', 'PaperSize', [20 10], 'PaperPosition', [0 0 20 10])
%print(fig6, 'Figure 6', '-dpdf');

fig7 = figure(7);
latqf1 = extractfield(Datosint, 'lat_QF');
latidx = find(latqf1 == 4);
gscatter(log(bbp1), log(bbpsat1),latqf1,'rgb','.',12);
title('lat qf');
xlabel('float'); ylabel('satellite');
set(fig7,'PaperUnits', 'centimeters', 'PaperSize', [20 10], 'PaperPosition', [0 0 20 10])
print(fig7, 'Figure 7', '-dpdf');
%% map 
filename = gunzip('gshhs_c.b.gz', tempdir);
shorelines = gshhs(filename{1});
delete(filename{1})
levels = [shorelines.Level];
land = (levels == 1);

figure,
hold on
worldmap([-90 -20],[0 360])
scatterm(lat1,lon1,[], bbp1,'k','filled'); colormap(jet)

geoshow(shorelines(land), 'FaceColor', [0.9 0.9 0.9])
hold off

%% [other]
figure,
worldmap([-90 -20],[0 360])
scatterm(latsat,lonsat,12,'k','filled'); colormap(jet)
geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])

figure,
worldmap([-90 90],[0 360])
scatterm(latsat,lonsat,50,'k','filled'); colormap(jet)
geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])

figure,
worldmap([-90 -20],[0 360])
scatterm(lat,lon,12,'k','filled'); colormap(jet)
geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
hold on
scatterm(latsat,lonsat,12,'r','filled'); colormap(jet)
b = find(lon>180);
lon(b) = lon(b) - 360;
figure,
scatter(lon,lonsat)
fields = fieldnames(Datoslid_SO)
bbpsat = extractfield(Datosint,'bbp');
bbp532sat = extractfield(Datosint,'bbp532');
bbp532sat = extractfield(Datosint,'bbp_532');
fields
bbp532sat = extractfield(Datosint,'bbp532');
fields = fieldnames(Datosint)
bbp532sat = extractfield(Datosint,'bbp532_sat_200');
bbpsat = extractfield(Datosint,'bbp532_sat_200');
bbp = extractfield(Datosint,'bbp_mld');
figure,scatter(bbp,bbpsat)
lsline
figure,hist(bbp)
figure,hist(bbpsat)
min(bbp)
min(bbpsat)
figure,scatter(log(bbp),log(bbpsat))
lsline
fitlm(bbp,bbpsat)
fitlm(log(bbp),log(bbpsat))

figure,
worldmap([-90 -20],[0 360])
scatterm(lat,lon,12,'k','filled'); colormap(jet)
geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
hold on
scatterm(latsat,lonsat,12,'r','filled'); colormap(jet)

figure,
worldmap([-90 90],[0 360])
scatterm(latsat,lonsat,50,'k','filled'); colormap(jet)
geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
%%
save Datosint2.mat


