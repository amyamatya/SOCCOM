% function [Datosint2] = float_lidar_matchups_II(Datoslid_SO,Datosint)
%cd('/Users/aamatya/Desktop/SOCCOM');
%load Datos_interpolated_12Mar2019.mat
%load Datos_lidar.mat

%% test
sat_times = extractfield(Datoslid_SO,'time_combined');
for xx = 1:10
    
    for ii = 1:length(Datosint2)
        display(ii)
        [a,b] = size(Datosint2(ii).time);
        for n = 1:b
            %display(n)
            float_month = month(Datosint2(ii).time(1,n));
            float_year = year(Datosint2(ii).time(1,n));
            float_day = day(Datosint2(ii).time(1,n));
            float_lat = Datosint2(ii).lat(1,n);
            float_lon = Datosint2(ii).lon(1,n);
            float_latQF = Datosint2(ii).lat_QF(1,n);
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
                    th = xx.*100;
                    outof_distance = find(thedistance > th);
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
            clearvars -except sat_lats_final sat_lons_final sat_bbp700_final thedays_final meand ii Datoslid_SO Datosint sat_times xx
        end
        %newidx = find(nonzeros(sat_bbp700_final));
        Datosint2(ii).res(xx).lat_sat = sat_lats_final;%(newidx);
        Datosint2(ii).res(xx).lon_sat = sat_lons_final;%(newidx);
        Datosint2(ii).res(xx).bbp700_sat = sat_bbp700_final;%(newidx);
        Datosint2(ii).res(xx).date_sat = thedays_final;%(newidx);
        Datosint2(ii).res(xx).meand_sat = meand;%(newidx);
        clearvars -except Datoslid_SO Datosint sat_times xx
    end
end
% Datosint2 = Datosint;
