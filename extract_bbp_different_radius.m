function [Datosint3] = extract_bbp_different_radius(Datosint2)
% Plots correlation between float and satellite bbp for varying radii
%% 
for ii = 1:length(Datosint2)
    Datosint3(ii).bbp1 = Datosint2(ii).res(1).bbp700_sat;
    Datosint3(ii).bbp2 = Datosint2(ii).res(2).bbp700_sat;
    Datosint3(ii).bbp3 = Datosint2(ii).res(3).bbp700_sat;
    Datosint3(ii).bbp4 = Datosint2(ii).res(4).bbp700_sat;
    Datosint3(ii).bbp5 = Datosint2(ii).res(5).bbp700_sat;
    Datosint3(ii).bbp6 = Datosint2(ii).res(6).bbp700_sat;
    Datosint3(ii).bbp7 = Datosint2(ii).res(7).bbp700_sat;
    Datosint3(ii).bbp8 = Datosint2(ii).res(8).bbp700_sat;
    Datosint3(ii).bbp9 = Datosint2(ii).res(9).bbp700_sat;
    Datosint3(ii).bbp10 = Datosint2(ii).res(10).bbp700_sat;
end

bbpsat1 = extractfield(Datosint3, 'bbp1');
bbpsat2 = extractfield(Datosint3, 'bbp2');
bbpsat3 = extractfield(Datosint3, 'bbp3');
bbpsat4 = extractfield(Datosint3, 'bbp4');
bbpsat5 = extractfield(Datosint3, 'bbp5');
bbpsat6 = extractfield(Datosint3, 'bbp6');
bbpsat7 = extractfield(Datosint3, 'bbp7');
bbpsat8 = extractfield(Datosint3, 'bbp8');
bbpsat9 = extractfield(Datosint3, 'bbp9');
bbpsat10 = extractfield(Datosint3, 'bbp10');
bbpfloat = extractfield(Datosint3, 'bbp_mld');
%%
diffdist_figure = figure;
hold on
sb1 = subplot(3,4,1); scatter(log(bbpfloat),log(bbpsat1),'.'); grid on; hold on; plot([-9 -4],[-8 -4],'k--');
set(gca, 'xlim',[-9 -4]); set(gca, 'ylim',[-8 -4]);
sb2 = subplot(3,4,2); scatter(log(bbpfloat),log(bbpsat2),'.'); hold on, plot([-9 -4],[-8 -4],'k--');grid on;
set(gca, 'xlim',[-9 -4]);set(gca, 'ylim',[-8 -4]);
sb3 = subplot(3,4,3);scatter(log(bbpfloat),log(bbpsat3),'.');hold on, plot([-9 -4],[-8 -4],'k--');grid on;
set(gca, 'xlim',[-9 -4]);set(gca, 'ylim',[-8 -4]);
sb4 = subplot(3,4,4);scatter(log(bbpfloat),log(bbpsat4),'.');hold on, plot([-9 -4],[-8 -4],'k--');grid on;
set(gca, 'xlim',[-9 -4]);set(gca, 'ylim',[-8 -4]);
sb5 = subplot(3,4,5);scatter(log(bbpfloat),log(bbpsat5),'.');hold on, plot([-9 -4],[-8 -4],'k--');grid on;
set(gca, 'xlim',[-9 -4]);set(gca, 'ylim',[-8 -4]);
sb6 = subplot(3,4,6);scatter(log(bbpfloat),log(bbpsat6),'.');hold on, plot([-9 -4],[-8 -4],'k--');grid on;
set(gca, 'xlim',[-9 -4]);set(gca, 'ylim',[-8 -4]);
sb7 = subplot(3,4,7);scatter(log(bbpfloat),log(bbpsat7),'.');hold on, plot([-9 -4],[-8 -4],'k--');grid on;
set(gca, 'xlim',[-9 -4]);set(gca, 'ylim',[-8 -4]);
sb8 = subplot(3,4,8);scatter(log(bbpfloat),log(bbpsat8),'.');hold on, plot([-9 -4],[-8 -4],'k--');grid on;
set(gca, 'xlim',[-9 -4]);set(gca, 'ylim',[-8 -4]);
sb9 = subplot(3,4,9);scatter(log(bbpfloat),log(bbpsat9),'.');hold on, plot([-9 -4],[-8 -4],'k--');grid on;
set(gca, 'xlim',[-9 -4]);set(gca, 'ylim',[-8 -4]);
sb10 = subplot(3,4,10);scatter(log(bbpfloat),log(bbpsat10),'.');hold on, plot([-9 -4],[-8 -4],'k--');grid on;
set(gca, 'xlim',[-9 -4]);set(gca, 'ylim',[-8 -4]);

addpath('/Users/aamatya/Desktop/suplabel');
sbl1 = suplabel('Float b_{bp}  (m^{-1}) ','x');
sbl2 = suplabel('Satellite b_{bp}  (m^{-1})','y');
sbl3 = suptitle('Satellite vs. Profiling Float Backscatter at 700nm');
set(sbl1, 'fontname','TimesNewRoman');
set(sbl2, 'fontname','TimesNewroman');
set(sbl3, 'fontname','timesnewroman');


lm1 = fitlm(log(bbpfloat),log(bbpsat1));
r2(1) = lm1.Rsquared(1).Ordinary;
rets(1) = lm1.NumObservations;
title(sb1, ['                  100km \newline \rm r^2 = ' num2str(r2(1)) ', ' num2str(rets(1)) ' retrievals']);

lm2 = fitlm(log(bbpfloat),log(bbpsat2));
r2(2) = lm2.Rsquared(1).Ordinary;
rets(2) = lm2.NumObservations;
title(sb2, ['                  200km \newline \rm r^2 = ' num2str(r2(2)) ', ' num2str(rets(2)) ' retrievals']);

lm3 = fitlm(log(bbpfloat),log(bbpsat3));
r2(3) = lm3.Rsquared(1).Ordinary;
rets(3) = lm3.NumObservations;
title(sb3, ['                  300km \newline \rm r^2 = ' num2str(r2(3)) ', ' num2str(rets(3)) ' retrievals']);

lm4 = fitlm(log(bbpfloat),log(bbpsat4));
r2(4) = lm4.Rsquared(1).Ordinary;
rets(4) = lm4.NumObservations;
title(sb4, ['                  400km \newline \rm r^2 = ' num2str(r2(4)) ', ' num2str(rets(4))  ' retrievals']);

lm5 = fitlm(log(bbpfloat),log(bbpsat5));
r2(5) = lm5.Rsquared(1).Ordinary;
rets(5) = lm5.NumObservations;
title(sb5, ['                  500km \newline \rm r^2 = ' num2str(r2(5)) ', ' num2str(rets(5)) ' retrievals']);

lm6 = fitlm(log(bbpfloat),log(bbpsat6));
r2(6) = lm6.Rsquared(1).Ordinary;
rets(6) = lm6.NumObservations;
title(sb6, ['                  600km \newline \rm r^2 = ' num2str(r2(6)) ', ' num2str(rets(6))  ' retrievals']);

lm7 = fitlm(log(bbpfloat),log(bbpsat7));
r2(7) = lm7.Rsquared(1).Ordinary;
rets(7) = lm7.NumObservations;
title(sb7, ['                  700km \newline \rm r^2 = ' num2str(r2(7)) ', ' num2str(rets(7)) ' retrievals']);

lm8 = fitlm(log(bbpfloat),log(bbpsat8));
r2(8) = lm8.Rsquared(1).Ordinary;
rets(8) = lm8.NumObservations;
title(sb8, ['                  800km \newline \rm r^2 = ' num2str(r2(8)) ', ' num2str(rets(8))  ' retrievals']);

lm9 = fitlm(log(bbpfloat),log(bbpsat9));
r2(9) = lm9.Rsquared(1).Ordinary;
rets(9) = lm9.NumObservations;
title(sb9, ['                  900km \newline \rm r^2 = ' num2str(r2(9)) ', ' num2str(rets(9))  ' retrievals']);

lm10 = fitlm(log(bbpfloat),log(bbpsat10));
r2(10) = lm10.Rsquared(1).Ordinary;
rets(10) = lm10.NumObservations;
title(sb10, ['                  1000km \newline \rm r^2 = ' num2str(r2(10)) ', ' num2str(rets(10)) ' retrievals']);

hold off
set(diffdist_figure,'PaperUnits', 'centimeters', 'PaperSize', [50 25], 'PaperPosition', [0 0 50 25])
print(diffdist_figure, 'diffdist', '-dpdf','-opengl');
%% r^2 plot
figure
plot(1:10, [r2(1) r2(2) r2(3) r2(4) r2(5) r2(6) r2(7) r2(8) r2(9) r2(10)]);
axis tight;
xlabel('distance');
ylabel('r^2');
title('dist vs. corrltn');





    
    
    
    
    