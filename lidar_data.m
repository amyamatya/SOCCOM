function [Datoslid_SO] = lidar_data

%% Datenum LIDAR Data
% (must re-add files to folder)
% 2012-2017 (Same years as float data) begins at row 135
cd('/Users/aamatya/Desktop/LIDAR');
addpath('/Users/aamatya/Desktop/natsortfiles');
%list = dir('/Users/aamatya/Desktop/LIDAR');    (original)
%sorted_list = natsortfiles({list(4:129).name});

%% 2012 - 2017
for n = 1:126
    cellfun(@load,sorted_list(n));
    Datoslid(n).name = sorted_list(n);
    Datoslid(n).time = profile_time;
    Datoslid(n).lat = lidar_latitude;
    Datoslid(n).lon = lidar_longitude;
    Datoslid(n).bbp = bbp2;
    Datoslid(n).bbp532 = bbp_532nm;
end
for n = 1:length(Datoslid)
    lat = Datoslid(n).lat;
    a = find(lat < -25);
    Datoslid_SO(n).name = Datoslid(n).name;
    Datoslid_SO(n).time = Datoslid(n).time(a);
    Datoslid_SO(n).lat = Datoslid(n).lat(a);
    Datoslid_SO(n).lon = Datoslid(n).lon(a);
    Datoslid_SO(n).bbp = Datoslid(n).bbp(a);
    Datoslid_SO(n).bbp532 = Datoslid(n).bbp532(a);
    clear a lat
end
%% 2006-2011
list = dir('/Users/aamatya/Desktop/LIDAR');
sorted_list = natsortfiles({list(4:137).name});

for n = 1:126
    Datoslid_SO(n+134) = Datoslid_SO(n);
end

for n = 1:length(sorted_list)
    cellfun(@load,sorted_list(n));
    Datoslid(n).name = sorted_list(n);
    Datoslid(n).time = profile_time;
    Datoslid(n).lat = lidar_latitude;
    Datoslid(n).lon = lidar_longitude;
    Datoslid(n).bbp = bbp2;
    Datoslid(n).bbp532 = bbp_532nm;
end

for n = 1:length(Datoslid)
    lat = Datoslid(n).lat;
    a = find(lat < -25);
    Datoslid_SO(n).name = Datoslid(n).name;
    Datoslid_SO(n).time = Datoslid(n).time(a);
    Datoslid_SO(n).lat = Datoslid(n).lat(a);
    Datoslid_SO(n).lon = Datoslid(n).lon(a);
    Datoslid_SO(n).bbp = Datoslid(n).bbp(a);
    Datoslid_SO(n).bbp532 = Datoslid(n).bbp532(a);
    clear a lat
end

for n = 1:length(Datoslid_SO)
    Datoslid_SO(n).year_simple = year(Datoslid_SO(n).time(1,:));
    Datoslid_SO(n).month_simple = month(Datoslid_SO(n).time(1,:));
    Datoslid_SO(n).time_combined = datenum(datetime(Datoslid_SO(n).year_simple,Datoslid_SO(n).month_simple,01));
end

clear sorted_list profile_time particle_depolarization n mes_sat list lidar_longitude...
    lidar_latitude kd_532nm Datoslid cross_polarization bbp_532nm bbp2 ans co_polarization
