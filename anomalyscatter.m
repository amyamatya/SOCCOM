%% Temp vs. Nitrate Anomaly Plot

trash1 = mapanoms_bbp;
trash2 = mapanoms;
trash3 = mapanoms_sst;
trash4 = maptimes;
trash5 = maptimes_sst

idx = find(isnan(trash2));
trash2(idx) = [];
trash3(idx) = [];
trash4(idx) = [];
trash5(idx) = [];

corrcoef(trash2, trash3)
p = polyfit(trash2, trash3, 1);
x1 = linspace(min(trash2), max(trash2), length(trash2));
y1 = polyval(p, x1);


fig19 = figure(19);
hold on
scatter(trash2, trash3,'.b');
plot(x1, y1,'r', 'linewidth',1.5)
axis tight
title('Nitrate V. Temperature Anomaly, r^2 = -0.8122');
xlabel('Nitrate');
ylabel('Temperature');
hold off

% print(fig19, '/Users/aamatya/Desktop/Summer 2019/19-tempnitscatter', '-dpdf');
%% Temp vs. Nitrate Scatter, by Month

for i = 1:12
    id = find(month(trash4) == i);
    trash(i).nit = trash2(id);
    trash(i).temp = trash3(id);
end

id1 = find(month(trash4) >= 10 | month(trash4) <= 2);
id2 = find(month(trash4) >= 3 & month(trash4) <= 9);
trash6nit = trash2(id1); trash6temp = trash3(id1);
trash7nit = trash2(id2); trash7temp = trash3(id2);

fig20 = figure(20);
hold on
scatter(trash7nit, trash7temp,40,  '.r');
scatter(trash6nit, trash6temp,40, '.b');
plot(x1, y1,'k', 'linewidth',1.5)
legend('Winter','Summer');
axis tight
title('Nitrate V. Temperature Anomaly, r^2 = -0.8122');
xlabel('Nitrate (\mum/kg)');
ylabel('Temperature (\circC)');
hold off
print(fig20, '/Users/aamatya/Desktop/Summer 2019/20-seasonscatter', '-dpdf');


% clear trash trash1 trash2 trash3 trash4 trash5 trash6nit trash6temp
% trash7nit trash7temp s1 s2