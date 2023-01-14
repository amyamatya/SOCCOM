%function [nitwoa18] = nitrate_test
% Extrac WOA18 nitrate data
cd('/Users/aamatya/Desktop/SOCCOM/nitdata');
list = dir('/Users/aamatya/Desktop/SOCCOM/nitdata');

for i = 4:length(list)
    display(i)
    nitrate_mean = ncread(list(i).name,'n_an');
    nitrate_mean = nitrate_mean(:,:,1:11);
    fillval_mean = ncreadatt(list(i).name,'n_an','_FillValue');
    nitrate_mean(nitrate_mean == fillval_mean) = NaN;
    nitrate_mean = nanmean(nitrate_mean,3);
    nitwoa18(i-3).nit = nitrate_mean;
end

lat = ncread(list(5).name,'lat');
lat = lat';
lat = repmat(lat,[360 1]);
nitwoa18(1).lat = lat;
lon = ncread(list(5).name,'lon');
lon = repmat(lon,[1 180]);
nitwoa18(1).lon = lon;


% Get data only south of 25S
a = find(lat<-25);
for n = 1:12
    nitwoa18(n).nit25S = nitwoa18(n).nit(a);
end
nitwoa18(1).lat25S = nitwoa18(1).lat(a);
nitwoa18(1).lon25S = nitwoa18(1).lon(a);

% Plots
filename = gunzip('gshhs_c.b.gz', tempdir);
shorelines = gshhs(filename{1});
delete(filename{1})
levels = [shorelines.Level];
land = (levels == 1);
% All data
figure,
for n = 1:12
    subplot(4,3,n)
    worldmap([-90 90],[0 360]),
    scatterm(lat(:),lon(:),[],nitwoa18(n).nit(:),'filled');
    title(num2str(n))
end
% South of 25S
lat25S = nitwoa18(1).lat25S;
lon25S = nitwoa18(1).lon25S;
figure,
for n = 1:12
    subplot(4,3,n)
    worldmap([-90 90],[0 360]),
    scatterm(lat25S(:),lon25S(:),[],nitwoa18(n).nit25S(:),'filled');
    title(num2str(n))
end

%%
for n = 1:12
    subplot(4,3,n)
    worldmap([-90 -25],[0 360]),
    scatterm(lat25S(:),lon25S(:),[],nitwoa18(n).nit25S(:),'filled');
    colormap(summer);colorbar;
    caxis([5 30]);
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
end