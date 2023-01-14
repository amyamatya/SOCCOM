cd('/Users/aamatya/Desktop/SOCCOM/sst_climatology');
list = dir('/Users/aamatya/Desktop/SOCCOM/sst_climatology');

%% Temperature Anomaly Map

for j = 4:length(list) % each month
    sst_mean = ncread(list(j).name,'t_an');
    sst_mean = sst_mean(:,1:65,1:11); %below 25S
    fillval_mean = ncreadatt(list(j).name,'t_an','_FillValue');
    sst_mean(sst_mean == fillval_mean) = NaN;
    sstdatas(1).overall(j-3).grid = nanmean(sst_mean,3); % Below x depth
end

for i = 1:length(Datosint2)
    Datosint2(i).time_combined = Datosint2(i).time(1,:);
    Datosint2(i).lat(Datosint2(i).lat_QF == 1) = NaN;
    Datosint2(i).lon(Datosint2(i).lat_QF == 1) = NaN;
    Datosint2(i).temperature_mld(Datosint2(i).lat_QF == 1) = NaN;
end

maplats_sst = [];
maplons_sst = [];
mapanoms_sst = [];
maptimes_sst = [];

for i = 1:length(Datosint2)
    for j = 1:length(Datosint2(i).lat)
        
        thelat = round(Datosint2(i).lat(j)) + 90;
        thelon = round(Datosint2(i).lon(j));
        thetime = Datosint2(i).month_simple(j);
  
        if isnan(thelat) | isnan(thelon)
            continue
        elseif thelon == 0
            thelon = 1;
        elseif thelon >= 360
            thelon = 360;
        end
        
        thebase = sstdatas(1).overall(thetime).grid(thelon, thelat);
         
        if isnan(thebase);
            continue
        end
        
        maplats_sst = [maplats_sst Datosint2(i).lat(j)];
        maplons_sst = [maplons_sst Datosint2(i).lon(j)];
        mapanoms_sst = [mapanoms_sst Datosint2(i).temperature_mld(j) - thebase];
        maptimes_sst = [maptimes_sst Datosint2(i).time_combined(j)];
    end
end
%% Temperature Anomaly Map
filename = gunzip('gshhs_c.b.gz', tempdir);
shorelines = gshhs(filename{1});
delete(filename{1})
levels = [shorelines.Level];
land = (levels == 1);

fig17 = figure(17)
worldmap([-90 -25],[0 360])
scatterm(maplats_sst ,maplons_sst,[], mapanoms_sst,'filled');colormap(redblue),colorbar
geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
title('Temperature Anomalies by float')
caxis([-7 7]);

% print(fig17, '/Users/aamatya/Desktop/Summer 2019/17-sstanomalymap','-dpdf');
%%

floatlat = extractfield(Datosint2, 'lat');
floatlon = extractfield(Datosint2, 'lon');
floattemp = extractfield(Datosint2, 'temperature_mld');

[newlon, idx] = sort(floatlon);
newlat = floatlat(idx);
newtemp = floattemp(idx);


plot(newlon, newtemp,'.')
xlabel('longitude');
ylabel('temp');



