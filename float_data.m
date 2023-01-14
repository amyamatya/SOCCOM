function [Datosint] = float_data
%% Data Extraction/Merging ([Datos]) - Merge float data in a single file to be used for NPP analyses

% Sort files by name order
cd('/Users/aamatya/Desktop/SOCCOM/float_data');
addpath('/Users/aamatya/Desktop/natsortfiles');
list=dir('/Users/aamatya/Desktop/SOCCOM/float_data');
sorted_list = natsortfiles({list(3:148).name}); 

% Loop through all files to extract desired variables
for n = 1:length(sorted_list)
    display(n)
    
    % Load data 
    cellfun(@load,sorted_list(n));
    
    % Combine date and time in a single variable
    mydate = datetime(datestr(FloatViz.mon_day_yr','mm/dd/yyyy'));
    mytime = timeofday(datetime(datestr(FloatViz.hh_mm')));
    mydatetime = mydate+mytime;
    time = datenum(mydatetime);
    [cols,rows]=size(FloatViz.Oxygen);
    time = repmat(time',cols,1);
    
    % Extract variables according to Quality Factor
    depth = FloatViz.Depth;
    nanes = find(FloatViz.Depth_QF>0);
    depth(nanes) = NaN;
    depth_QF = FloatViz.Depth_QF;
    clear nanes
    
    oxygen = FloatViz.Oxygen;
    nanes = find(FloatViz.Oxygen_QF>0);
    oxygen(nanes) = NaN;
    oxygen_QF = FloatViz.Oxygen_QF;
    oxygen  = smoothdata(oxygen,1,'movmedian',7);
    clear nanes
    
    lat = FloatViz.Lat;
    nanes = find(FloatViz.Lat_QF>4);
    lat(nanes) = NaN;
    lat_QF = FloatViz.Lat_QF;
    lon = FloatViz.Lon;
    lon(nanes) = NaN;
    clear nanes
    
    nitrate = FloatViz.Nitrate;
    nanes = find(FloatViz.Nitrate_QF>0);
    nitrate(nanes) = NaN;
    nitrate_QF = FloatViz.Nitrate_QF;
    nitrate  = smoothdata(nitrate,1,'movmedian',7);
    clear nanes
    
    poc = FloatViz.POC;
    nanes = find(FloatViz.POC_QF>1);
    poc(nanes) = NaN;
    poc_QF = FloatViz.POC_QF;
    poc(poc == -1.0000e+10) = NaN;
    poc  = smoothdata(poc,1,'movmedian',7);
    clear nanes
    
    bbp = FloatViz.b_bp700;
    nanes = find(FloatViz.b_bp700_QF>1);
    bbp(nanes) = NaN;
    bbp(bbp==-1.0000e+10)=NaN;
    bbp  = smoothdata(bbp,1,'movmedian',7);
    bbp_QF = FloatViz.b_bp700_QF;
    clear nanes
    
    chl = FloatViz.Chl_a_corr;
    nanes = find(FloatViz.Chl_a_corr_QF>1);
    chl(nanes) = NaN;
    chl(chl == -1.0000e+10) = NaN;
    chl_QF = FloatViz.Chl_a_corr_QF;
    chl  = smoothdata(chl,1,'movmedian',7);
    clear nanes
    
    chl_uncor = FloatViz.Chl_a;
    chl_uncor(chl_uncor == -1.0000e+10) = NaN;
    nanes = find(FloatViz.Chl_a_QF>1);
    chl_uncor(nanes) = NaN;
    chl_uncor_QF = FloatViz.Chl_a_QF;
    chl_uncor  = smoothdata(chl_uncor,1,'movmedian',7);
    clear nanes
    
    temperature = FloatViz.Temperature;
    nanes = find(FloatViz.Temperature_QF>0);
    temperature(nanes) = NaN;
    temperature = smoothdata(temperature,1,'movmedian',7);
    temperature_QF = FloatViz.Temperature_QF;
    clear nanes
    
    pressure = FloatViz.Pressure;
    nanes = find(FloatViz.Pressure_QF>0);
    pressure(nanes) = NaN;
    pressure = smoothdata(pressure,1,'movmedian',7);
    pressure_QF = FloatViz.Pressure_QF;
    clear nanes
    
    salinity = FloatViz.Salinity;
    nanes = find(FloatViz.Salinity_QF>0);
    salinity(nanes) = NaN;
    salinity = smoothdata(salinity,1,'movmedian',7);
    salinity_QF = FloatViz.Salinity_QF;
    clear nanes
    
    sigma_theta = FloatViz.Sigma_theta;
    nanes = find(FloatViz.Sigma_theta_QF>0);
    sigma_theta(nanes) = NaN;
    sigma_theta  = smoothdata(sigma_theta,1,'movmedian',7);
    sigma_theta_QF = FloatViz.Sigma_theta_QF;
    clear nanes
    
    % Carbonate system variables
    varnames = fieldnames(FloatViz);
    if ~isempty(find(strcmp('pHinsitu',varnames) == 1)) % Check if float has pHinsitu
        
        pHinsitu = FloatViz.pHinsitu;
        nanes = find(FloatViz.pHinsitu_QF>0);
        pHinsitu  = smoothdata(pHinsitu,1,'movmedian',7);
        pHinsitu(nanes) = NaN;
        pHinsitu_QF = FloatViz.pHinsitu_QF;
        clear nanes
        
        pH25C = FloatViz.pH25C;
        nanes = find(FloatViz.pH25C_QF>0);
        pH25C = smoothdata(pH25C,1,'movmedian',7);
        pH25C(nanes) = NaN;
        pH25C_QF = FloatViz.pH25C_QF;
        clear nanes
        
        TALK_LIAR = FloatViz.TALK_LIAR;
        nanes = find(FloatViz.TALK_LIAR_QF>0);
        TALK_LIAR = smoothdata(TALK_LIAR,1,'movmedian',7);
        TALK_LIAR(nanes) = NaN;
        TALK_LIAR_QF = FloatViz.TALK_LIAR_QF;
        clear nanes
        
        DIC_LIAR = FloatViz.DIC_LIAR;
        nanes = find(FloatViz.DIC_LIAR_QF>0);
        DIC_LIAR = smoothdata(DIC_LIAR,1,'movmedian',7);
        DIC_LIAR(nanes) = NaN;
        DIC_LIAR_QF = FloatViz.DIC_LIAR_QF;
        clear nanes
        
        pCO2_LIAR = FloatViz.pCO2_LIAR;
        nanes = find(FloatViz.pCO2_LIAR_QF>0);
        pCO2_LIAR = smoothdata(pCO2_LIAR,1,'movmedian',7);
        pCO2_LIAR(nanes) = NaN;
        pC02_LIAR_QF = FloatViz.pCO2_LIAR;
        clear nanes
    else
    end
    
    % Get years of analysis for current float in the loop
    fname = strsplit(string(sorted_list(n)),'.');
    final_name = char(strcat('f',fname(1)));
    
    Data.(final_name).fname = final_name;
    Data.(final_name).time = time;
    Data.(final_name).depth = depth;
    Data.(final_name).depth_QF = depth_QF;
    Data.(final_name).lat = lat;
    Data.(final_name).lat_QF = lat_QF;
    Data.(final_name).lon = lon;
    Data.(final_name).oxygen = oxygen;
    Data.(final_name).oxygen_QF = oxygen_QF;
    Data.(final_name).chl = chl;
    Data.(final_name).chl_QF = chl_QF;
    Data.(final_name).chl_uncor = chl_uncor;
    Data.(final_name).chl_uncor_QF = chl_uncor_QF;
    Data.(final_name).salinity = salinity;
    Data.(final_name).salinity_QF = salinity_QF;
    Data.(final_name).temperature = temperature;
    Data.(final_name).temperature_QF = temperature_QF;
    Data.(final_name).pressure = pressure;
    Data.(final_name).pressure_QF = pressure_QF;
    Data.(final_name).sigma_theta = sigma_theta;
    Data.(final_name).sigma_theta_QF = sigma_theta_QF;
    Data.(final_name).nitrate = nitrate;
    Data.(final_name).nitrate_QF = nitrate_QF;
    Data.(final_name).poc = poc;
    Data.(final_name).poc_QF = poc_QF;
    Data.(final_name).bbp = bbp;
    Data.(final_name).bbp_QF = bbp_QF;
    
    if isempty(find(strcmp('pHinsitu',varnames) == 1)) == 0
        Data.(final_name).pHinsitu = pHinsitu;
        Data.(final_name).pH25C = pH25C;
        Data.(final_name).TALK_LIAR = TALK_LIAR;
        Data.(final_name).DIC_LIAR = DIC_LIAR;
        Data.(final_name).pCO2_LIAR = pCO2_LIAR;
    else
        Data.(final_name).pHinsitu = NaN;
        Data.(final_name).pH25C = NaN;
        Data.(final_name).TALK_LIAR = NaN;
        Data.(final_name).DIC_LIAR = NaN;
        Data.(final_name).pCO2_LIAR = NaN;
    end
end

Datos = struct2array(Data);

for n = 1:length(Datos)
    Datos(n).numdays_sampled = Datos(n).time(1,end) - Datos(n).time(1,1);
end

%% Interpolation ([Datosint]) - Interpolate float data to 2000m with a 1m vertical resolution

%load('npp_datos.mat')
for i = 1:length(Datos)
    display(i) 
    xq = [1:2000]';
    [a b] = size(Datos(i).chl);
    for n = 1:b
        x = Datos(i).depth(:,n);
        a = find(isnan(x));
        x(a) = [];
        c = Datos(i).chl(:,n);
        c(a) = [];
        cu = Datos(i).chl_uncor(:,n);
        cu(a) = [];
        p = Datos(i).poc(:,n);
        p(a) = [];
        t = Datos(i).temperature(:,n);
        t(a) = [];
        s = Datos(i).salinity(:,n);
        s(a) = [];
        d = Datos(i).sigma_theta(:,n);
        d(a) = [];
        nit = Datos(i).nitrate(:,n);
        nit(a) = [];
        ox = Datos(i).oxygen(:,n);
        ox(a) = [];
        bbp = Datos(i).bbp(:,n);
        bbp(a) = [];
        if length(x) == length(unique(x)) 
            cq(:,n) = interp1(x,c,xq);
            cuq(:,n) = interp1(x,cu,xq);
            pq(:,n) = interp1(x,p,xq);
            tq(:,n) = interp1(x,t,xq);
            sq(:,n) = interp1(x,s,xq);
            dq(:,n) = interp1(x,d,xq);
            nitq(:,n) = interp1(x,nit,xq);
            oxq(:,n) = interp1(x,ox,xq);
            bbpq(:,n) = interp1(x,bbp,xq);
        else
            pre_cq(1:2000,1) = NaN;    
            cq(:,n) = pre_cq;
            pre_cuq(1:2000,1) = NaN;
            cuq(:,n) = pre_cuq;
            pre_pq(1:2000,1) = NaN;
            pq(:,n) = pre_pq;
            pre_tq(1:2000,1) = NaN;
            tq(:,n) = pre_tq;
            pre_sq(1:2000,1) = NaN;
            sq(:,n) = pre_sq;
            pre_dq(1:2000,1) = NaN;
            dq(:,n) = pre_dq;
            pre_nitq(1:2000,1) = NaN;
            nitq(:,n) = pre_nitq;
            pre_oxq(1:2000,1) = NaN;
            oxq(:,n) = pre_oxq;
            pre_bbpq(1:2000,1) = NaN;
            bbpq(:,n) = pre_bbpq;
        end
    end
    Datosint(i).fname = Datos(i).fname;
    Datosint(i).time = Datos(i).time;
    Datosint(i).lat = Datos(i).lat;
    Datosint(i).lat_QF = Datos(i).lat_QF;
    Datosint(i).lon = Datos(i).lon;
    Datosint(i).chl = cq;
    Datosint(i).chl_uncor = cuq;
    Datosint(i).poc = pq;
    Datosint(i).temperature = tq;
    Datosint(i).salinity = sq;
    Datosint(i).density = dq;
    Datosint(i).nitrate = nitq;
    Datosint(i).oxygen = oxq;
    Datosint(i).bbp = bbpq;
    Datosint(i).numdays_sampled = Datos(i).numdays_sampled;
    %Datosint(i).par = Datos(i).par;
    %Datosint(i).mld = Datos(i).mld;
    clearvars -except Datos Datosint
end

%for qwe = 1: length(Datosint)
    
%save Datosint.mat;
%SOCCOM_time = extractfield(Datosint, 'time');
%SOCCOM_bbp = extractfield(Datosint, 'bbp');
%save ('SOCCOMfinal.mat', 'SOCCOM_time','SOCCOM_bbp');
%% Map Plot - Plotting ancillary files

lat = Datosint(1).lat % latitudinal position of the data
lon = Datosint(1).lon; % longitudinal position of the data
x = Datosint(1).chl; % the variable that you want to plot (use a color code (e.g., 'k') 
%if you just want to see the position of the data)
for i = 1:length(Datosint)
    [a b] = size(Datosint(i).chl);
    for n = 1:b
        chl = Datosint(i).chl(1:20,n);
        chl_20m(n) = nanmean(chl);
    end
    Datosint(i).chl_20m = chl_20m;
    clear chl chl_20m
end
%figure,plot(Datos(1).depth(1:20,1),'o')
%figure
%worldmap([-90 -20],[0 360])
%scatterm(lats,lons,[],chl_mlds,'filled'); colormap(jet)
%geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
%% MLD Algorithm - input insitu temp, salinity, depths/pressures ... output MLD, isothermal layer depth, sorted profile of potential density
addpath('/Users/aamatya/Desktop/SOCCOM/SW_Functions');
cd('/Users/aamatya/Desktop/SOCCOM');
[c d] = size(Datos);
for i = 1:d
    [a b] = size(Datos(i).temperature);
    for j = 1:b
        [mld_out_Datos,ild_out_Datos,sig_theta_Datos] = mld_dbm(Datos(i).temperature(:,j),Datos(i).salinity(:,j),Datos(i).depth(:,j),0);
        mld(j) = mld_out_Datos;
    end
    Datos(i).mld = mld;
    Datosint(i).mld = mld;
    clear mld mld_out
end
%% Averaging MLD for Chl
for i = 1:length(Datosint)
    [a b] = size(Datosint(i).chl);
    for n = 1:b
        rmld = round(Datosint(i).mld(:,n));
        if isnan(rmld) == 0
            chl = Datosint(i).chl(1:rmld,n);
            chl_mld(n) = nanmean(chl);
        else
            chl_mld(n) = NaN;
        end
    end
    Datosint(i).chl_mld = chl_mld;
    clear chl chl_mld mld
end

lat_QFs = extractfield(Datosint, 'lat_QF');
lats = extractfield(Datosint, 'lat');
lons = extractfield(Datosint, 'lon');
chl_mlds = extractfield(Datosint, 'chl_mld');
f = find(lat_QFs == 1);
lats(f) = NaN;
lons(f) = NaN;
chl_mlds(f) = NaN;

%lat = Datosint(1).lat % latitudinal position of the data
%lon = Datosint(1).lon; % longitudinal position of the data
%x = Datosint(1).chl_mld; %
%figure
%worldmap([-90 -20],[0 360])
%scatterm(lats,lons,[], chl_mlds,'filled'); colormap(jet)
%geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
%colorbar;
%% MLD for Chl (uncorrected)
for i = 1:length(Datosint)
    [a b] = size(Datosint(i).chl_uncor);
    for n = 1:b
        rmld = round(Datosint(i).mld(:,n));
        if isnan(rmld) == 0
            chl_uncor = Datosint(i).chl_uncor(1:rmld,n);
            chl_uncor_mld(n) = nanmean(chl_uncor);
        else
            chl_uncor_mld(n) = NaN;
        end
    end
    Datosint(i).chl_uncor_mld = chl_uncor_mld;
    clear chl_uncor chl_uncor_mld mld
end

%figure
%worldmap([-90 -20],[0 360])
%scatterm(lats,lons,[],chl_uncor_mlds,'filled'); colormap(jet)
%geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
%colorbar;
%title('Chlorophyl Uncorrected');
%% MLD for POC
for i = 1:length(Datosint)
    [a b] = size(Datosint(i).poc);
    for n = 1:b
        rmld = round(Datosint(i).mld(:,n));
        if isnan(rmld) == 0
            poc = Datosint(i).poc(1:rmld,n);
            poc_mld(n) = nanmean(poc);
        else
            poc_mld(n) = NaN;
        end
    end
    Datosint(i).poc_mld = poc_mld;
    clear poc poc_mld mld
end

poc_mlds = extractfield(Datosint, 'poc_mld');
poc_mlds(f) = NaN;

%figure
%worldmap([-90 -20],[0 360])
%scatterm(lats,lons,[], poc_mlds,'filled'); colormap(jet)
%geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
%colorbar;
%title('POC');
%% MLD for Temp

for i = 1:length(Datosint)
    [a b] = size(Datosint(i).temperature);
    for n = 1:b
        rmld = round(Datosint(i).mld(:,n));
        if ~isnan(rmld)
            temperature = Datosint(i).temperature(1:rmld,n);
            temperature_mld(n) = nanmean(temperature);
        else
            temperature_mld(n) = NaN;
        end
    end
    Datosint(i).temperature_mld = temperature_mld;
    clear temperature temperature_mld mld
end

temperature_mlds = extractfield(Datosint, 'temperature_mld');
temperature_mlds(f) = NaN;

%figure
%worldmap([-90 -20],[0 360])
%scatterm(lats,lons,[], temperature_mlds,'filled'); colormap(jet)
%geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
%colorbar;
%title('Temperature');
%% MLD for Salinity
for i = 1:length(Datosint)
    [a b] = size(Datosint(i).salinity);
    for n = 1:b
        rmld = round(Datosint(i).mld(:,n));
        if isnan(rmld) == 0
            salinity = Datosint(i).salinity(1:rmld,n);
            salinity_mld(n) = nanmean(salinity);
        else
            salinity_mld(n) = NaN;
        end
    end
    Datosint(i).salinity_mld = salinity_mld;
    clear salinity salinity_mld mld
end

salinity_mlds = extractfield(Datosint, 'salinity_mld');
salinity_mlds(f) = NaN;

%figure
%worldmap([-90 -20],[0 360])
%scatterm(lats,lons,[], salinity_mlds,'filled'); colormap(jet)
%geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
%colorbar;
%title('Salinity');
%% MLD for Density
for i = 1:length(Datosint)
    [a b] = size(Datosint(i).density);
    for n = 1:b
        rmld = round(Datosint(i).mld(:,n));
        if isnan(rmld) == 0
            density = Datosint(i).density(1:rmld,n);
            density_mld(n) = nanmean(density);
        else
            density_mld(n) = NaN;
        end
    end
    Datosint(i).density_mld = density_mld;
    clear density density_mld mld
end

density_mlds = extractfield(Datosint, 'density_mld');
density_mlds(f) = NaN;

%figure
%worldmap([-90 -20],[0 360])
%scatterm(lats,lons,[], density_mlds,'filled'); colormap(jet)
%geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
%colorbar;
%title('Density');
%% MLD for Nitrate
for i = 1:length(Datosint)
    [a b] = size(Datosint(i).nitrate);
    for n = 1:b
        rmld = round(Datosint(i).mld(:,n));
        if isnan(rmld) == 0
            nitrate = Datosint(i).nitrate(1:rmld,n);
            nitrate_mld(n) = nanmean(nitrate);
        else
            nitrate_mld(n) = NaN;
        end
    end
    Datosint(i).nitrate_mld = nitrate_mld;
    clear nitrate nitrate_mld mld
end

nitrate_mlds = extractfield(Datosint, 'nitrate_mld');
nitrate_mlds(f) = NaN;

%figure
%worldmap([-90 -20],[0 360])
%scatterm(lats,lons,[], nitrate_mlds,'filled'); colormap(jet)
%geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
%colorbar;
%title('Nitrate');
%% Mld for Oxygen
for i = 1:length(Datosint)
    [a b] = size(Datosint(i).oxygen);
    for n = 1:b
        rmld = round(Datosint(i).mld(:,n));
        if isnan(rmld) == 0
            oxygen = Datosint(i).oxygen(1:rmld,n);
            oxygen_mld(n) = nanmean(oxygen);
        else
            oxygen_mld(n) = NaN;
        end
    end
    Datosint(i).oxygen_mld = oxygen_mld;
    clear oxygen oxygen_mld mld
end
oxygen_mlds = extractfield(Datosint, 'oxygen_mld');
oxygen_mlds(f) = NaN;

%figure
%worldmap([-90 -20],[0 360])
%scatterm(lats,lons,[], oxygen_mlds,'filled'); colormap(jet)
%geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
%colorbar;
%title('Oxygen');
%% MLD for BBP
for i = 1:length(Datosint)
    [a b] = size(Datosint(i).bbp);
    for n = 1:b
        rmld = round(Datosint(i).mld(:,n));
        if isnan(rmld) == 0
            bbp = Datosint(i).bbp(1:rmld,n);
            bbp_mld(n) = nanmean(bbp);
        else
            bbp_mld(n) = NaN;
        end
    end
    Datosint(i).bbp_mld = bbp_mld;
    clear bbp bbp_mld mld
end


bbp_mlds = extractfield(Datosint, 'bbp_mld');
f = find(lat_QFs >= 1);
bbp_mlds(f) = NaN;

%figure
%worldmap([-90 -20],[0 360])
%scatterm(lats,lons, [], bbp_mlds,'filled'); colormap(jet)
%geoshow(shorelines(land),  'FaceColor', [0.9 0.9 0.9])
%title('BBP');

for a = 1:length(Datosint)
    Datosint(a).year = year(Datosint(a).time);
    Datosint(a).month = month(Datosint(a).time);
    Datosint(a).day = day(Datosint(a).time);
end

