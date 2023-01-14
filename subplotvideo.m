%% Float + Satellite all, monthly (see clim_test)
titles = {'January','February','March','April','May','June','July','August','September','October','November','December'};
for n = 1:12
    fig(n) = figure(n+19);
    hold on
    worldmap([-90 -20],[0 360]);
    trash1 = scatterm(clim(n).lat_sat,clim(n).lon_sat,0.5, [1 0 0],'filled');
    trash2 = scatterm(clim(n).lat,clim(n).lon,3,[0 0 1], 'filled');
    
    [~, hObj] = legend([trash1 trash2],'Satellite','Float', 'location','northeastoutside');
    hL=findobj(hObj,'type','line');
    set(hL,'linewidth',5,'location','northwest');
    
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9]);
    title(titles(n), 'Position',[0 8.3e6]);
    hold off
    f(n) = getframe(fig(n));
end
close all

for i = 1:12
    f(i+12) = f(i);
    f(i+24) = f(i);
    f(i+36) = f(i);
    f(i+48) = f(i);
    f(i+60) = f(i);
    f(i+72) = f(i);
    f(i+84) = f(i);
    f(i+96) = f(i);
end

v = VideoWriter('anomvideo.mp4', 'MPEG-4');
v.FrameRate = 4;
v.Quality = 100;
open(v)
for i = 1:96
    writeVideo(v, f(i));
end
close(v)
clear f;
%% monthly satellite bbp anomaly (see bbpanomaly)

for k = 1:12
    fig(k) = figure(k);
    worldmap([-90 -25],[0 360])
    scatterm(yearlybbp(k).satlat, yearlybbp(k).satlon, [], yearlybbp(k).satanom,'filled'); colormap(redblue),colorbar
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
    caxis([-0.0006 0.0006]);
    title(titles(k), 'Position',[0 8e6]);
    f(k) = getframe(fig(k));
end
close all

for i = 1:12
    f(i+12) = f(i);
    f(i+24) = f(i);
    f(i+36) = f(i);
    f(i+48) = f(i);
    f(i+60) = f(i);
    f(i+72) = f(i);
    f(i+84) = f(i);
    f(i+96) = f(i);
end
v = VideoWriter('satanomvideo.mp4', 'MPEG-4');
v.FrameRate = 4;
v.Quality = 100;
open(v)
for i = 1:96
    writeVideo(v, f(i));
end
close(v)
clear f;
%% monthly satellite bbp total (see bbpanomaly)
% figure,
for k = 1:12
    fig(k) = figure(k);
    worldmap([-90 -25],[0 360])
    scatterm(yearlybbp(k).satlat_tot, yearlybbp(k).satlon_tot, [], yearlybbp(k).satanom_tot,'filled'); colormap(summer),colorbar
    geoshow(shorelines(land),'FaceColor', [0.9 0.9 0.9])
    title(titles(k), 'Position',[0 8e6]);
    %caxis([0 .004]);
%     f(k) = getframe(fig(k));
end
close all;
%%
for i = 1:12
    f(i+12) = f(i);
    f(i+24) = f(i);
    f(i+36) = f(i);
    f(i+48) = f(i);
    f(i+60) = f(i);
    f(i+72) = f(i);
    f(i+84) = f(i);
    f(i+96) = f(i);
end

v = VideoWriter('sattotal.mp4', 'MPEG-4');
v.FrameRate = 4;
v.Quality = 100;
open(v)
for i = 1:96
    writeVideo(v, f(i));
end
close(v)
clear f;
%% nitrate anom + bbp total
addpath('/Users/aamatya/Desktop/suplabel');

for i = 1%:12
    fig(i) = figure(i);
%     set(gca, 'Position',[100 100 800 300]);
    trash1 = subplot(1,2,1)
    worldmap([-90 -25],[0 360])
    scatterm(nitdata(4).monthlymap(i).lat ,nitdata(4).monthlymap(i).lon,[], nitdata(4).monthlymap(i).anom,'filled'); colormap(trash1, redblue), trash = colorbar;
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
    axis equal;
    caxis([-10 10])
    
    trash2 = subplot(1,2,2)
    worldmap([-90 -25],[0 360])
    scatterm(yearlybbp(i).satlat_tot, yearlybbp(i).satlon_tot, [], yearlybbp(i).satanom_tot,'filled'); colormap(trash2, summer) ,colorbar;
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
    axis equal;
    caxis([0 .004]);
    
    suptitle(titles(i));
%     f(i) = getframe(fig(i));
end

close all;

for i = 1:12
    f(i+12) = f(i);
    f(i+24) = f(i);
    f(i+36) = f(i);
    f(i+48) = f(i);
    f(i+60) = f(i);
    f(i+72) = f(i);
    f(i+84) = f(i);
    f(i+96) = f(i);
end

v = VideoWriter('satcompare.mp4', 'MPEG-4');
v.FrameRate = 2;
v.Quality = 100;
open(v)
for i = 1:96
    writeVideo(v, f(i));
end
close(v)
clear f;

%% nitrate total + bbp total

for i = 1:12
    fig(i) = figure(i);
%     set(gca, 'Position',[100 100 800 300]);
    trash1 = subplot(1,2,1)
    worldmap([-90 -25],[0 360]),
    scatterm(lat25S(:),lon25S(:),[],nitwoa18(n).nit25S(:),'filled');
    colormap(trash1, summer);colorbar;
    caxis([5 30]);
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
    axis equal;
    
    trash2 = subplot(1,2,2)
    worldmap([-90 -25],[0 360])
    scatterm(yearlybbp(i).satlat_tot, yearlybbp(i).satlon_tot, [], yearlybbp(i).satanom_tot,'filled'); colormap(trash2, summer); colorbar
    geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
    axis equal;
    caxis([0 .004]);
    
    suptitle(titles(i));
    f(i) = getframe(fig(i));
end

close all;

for i = 1:12
    f(i+12) = f(i);
    f(i+24) = f(i);
    f(i+36) = f(i);
    f(i+48) = f(i);
    f(i+60) = f(i);
    f(i+72) = f(i);
    f(i+84) = f(i);
    f(i+96) = f(i);
end

v = VideoWriter('totalcompare.mp4', 'MPEG-4');
v.FrameRate = 2;
v.Quality = 100;
open(v)
for i = 1:96
    writeVideo(v, f(i));
end
close(v)
clear f;
