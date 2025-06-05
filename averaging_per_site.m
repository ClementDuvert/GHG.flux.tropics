%% Function to average out values from sites with several measurements (using either means or medians).
% For large rivers, this function merges nearby sites within 5km.
% For lakes/reservoirs, this function averages out concentration and flux values for each system.
%
function[T,Ts]=averaging_per_site(T,choice)
%
% Inputs
%  - T must be a table with all variables listed as per txt file, and must have gone through previous functions such as 'unit_conversion' and 'climate_classification'
%  - choice is 'streams/rivers' or 'lakes/reservoirs'
%  - requires the following files in the input data folder (for lakes):
% 'lake_HydroLAKES_2025.txt' (HydroLAKES data)
% 'lake_wind_2025.csv'       (WorldClim2 wind data for each system)
% 'lake_manual_2025.csv'     (surface area estimates for systems not in HydroLAKES)
%
% Outputs
%  - T is the unchanged input table
%  - Ts is the reduced table with one row per site and mean/median values

disp(' '); disp('starting calculation of means/medians per site...'); disp(' ')


%% 0. Input datasets + initial calculations for lakes

%constants
CO2_sat=415; %umol/mol or ppmv
CH4_sat=1.87; %umol/mol or ppmv
N2O_sat=0.332; %umol/mol or ppmv  in 2020
gas_constant=8.3144626181; %m3 Pa mol-1 K-1
threshold_km=5;  %threshold for large river location merge
SO_large=7; %'large river' stream order cutoff
minpts=2; %for merging

if isequal(choice,'lakes/reservoirs')

    %import dataset with HydroLAKES data and add area, depth etc. to table
    SA=readtable('lake_HydroLAKES_2025.txt');
    [isMatch,idxInSA]=ismember(T.Lake_idx,SA.group_idx);
    T.hylak_id=NaN(height(T),1);
    T.hylak_lake_name=strings(height(T),1);
    T.hylak_lake_area=NaN(height(T),1);
    T.hylak_vol_total=NaN(height(T),1);
    T.hylak_depth_avg=NaN(height(T),1);
    T.hylak_id(isMatch)=SA.hylak_id(idxInSA(isMatch));
    T.hylak_lake_name(isMatch)=SA.hylak_lake_name(idxInSA(isMatch));
    T.hylak_lake_area(isMatch)=SA.hylak_lake_area(idxInSA(isMatch));
    T.hylak_vol_total(isMatch)=SA.hylak_vol_total(idxInSA(isMatch));
    T.hylak_depth_avg(isMatch)=SA.hylak_depth_avg(idxInSA(isMatch));

    % import dataset for lakes that are not in HydroLAKES and add area, depth etc. to table
    SB=readtable('lake_manual_2025.csv');
    [isMatchB,idxInSB]=ismember(T.Lake_idx,SB.group_idx_manual);
    T.lake_area_manual=NaN(height(T),1); %km2
    T.vol_total_manual=NaN(height(T),1); %Mm3
    T.depth_avg_manual=NaN(height(T),1); %m
    T.lake_area_manual(isMatchB)=SB.lake_area_km2(idxInSB(isMatchB));
    T.vol_total_manual(isMatchB)=SB.vol_total_Mm3(idxInSB(isMatchB));
    T.depth_avg_manual(isMatchB)=SB.depth_avg_m(idxInSB(isMatchB));

    % consolidate surface areas (in km2)
    T.area_consolidated=T.hylak_lake_area; %preferred area data = HydroLAKES (km2)
    for s=1:size(T,1)
        if isnan(T.area_consolidated(s))
            T.area_consolidated(s)=T.SurfaceArea_m2(s)/1E6; %if no HydroLAKES we take the area reported during extraction (m2 to km2)
            if isnan(T.SurfaceArea_m2(s))
                T.area_consolidated(s)=T.lake_area_manual(s); %if no area reported we use manual delineation data
            end
        end
    end

    %import wind speed data and incorporate into table
    SC=readtable('lake_wind_2025.csv');
    [isMatchC,idxInSC]=ismember(T.Lake_idx,SC.lake_ID);
    T.wind=NaN(height(T),1);
    T.wind(isMatchC)=SC.wind_worldclim(idxInSC(isMatchC));

    % calculate Schmidt numbers
    T.Sc_CO2=1923.6-125.06.*T.Temp_gapfill+4.3773.*T.Temp_gapfill.^2-0.085681.*T.Temp_gapfill.^3+0.00070284.*T.Temp_gapfill.^4; %Jahne et al (1987) / Wanninkhof (2014)
    T.Sc_CH4=1909.4-120.78.*T.Temp_gapfill+4.1555.*T.Temp_gapfill.^2-0.080578.*T.Temp_gapfill.^3+0.00065777.*T.Temp_gapfill.^4; %Jahne et al (1987) / Wanninkhof (2014)
    T.Sc_N2O=2141.2-152.56.*T.Temp_gapfill+5.8963.*T.Temp_gapfill.^2-0.12411.*T.Temp_gapfill.^3+0.0010655.*T.Temp_gapfill.^4; %Jahne et al (1987) / Wanninkhof (2014)

    % calculate k600 based on Wanninkhof 1992; Cole and Caraco 1998; Vachon and Prairie 2013
    for a=1:size(T,1)
        w1=0.45*T.wind(a)^1.64; %in cm/h
        T.k600_modelled_W92(a)=w1/100*24; %in m/d
        w2=2.07+0.215*T.wind(a)^1.7; %in cm/h
        T.k600_modelled_CC98(a)=w2/100*24; %in m/d
        w3=2.51+(1.48*T.wind(a))+(0.39*T.wind(a)*log10(T.area_consolidated(a))); %in cm/h
        T.k600_modelled_VP13(a)=w3/100*24; %in m/d
        clear w1 w2 w3
    end

    %consolidate k600 - we keep the mean of CC98 and VP13 if no reported value (see flux_calculation_lakes.m)
    T.k600_CC98VP13=mean([T.k600_modelled_CC98,T.k600_modelled_VP13],2);
    T.k600_consolidated=T.k600_converted; %we always prefer reported k600
    T.k600_consolidated(isnan(T.k600_converted))=T.k600_CC98VP13(isnan(T.k600_converted)); %if no reported k600, use modelled one

    % calculate kCO2, kCH4 and kN2O
    gases={'CO2','CH4','N2O'};
    models={'W92','CC98','VP13'};
    for i=1:length(models)
        for j=1:length(gases)
            gas=gases{j};
            model=models{i};
            T.(['k' lower(gas) '_' model])=T.(['k600_modelled_' model]).*(T.(['Sc_' gas])/600).^-0.5;
        end
    end

    % calculate Henry constants and atmospheric concentrations
    T.Temp_Kelvin=273.15+T.Temp_gapfill;
    T.Kh_CO2=10.^-(-(108.3865+0.01985076.*(T.Temp_Kelvin)-6919.53./(T.Temp_Kelvin)-40.4515.*log10(T.Temp_Kelvin)+669365./(T.Temp_Kelvin).^2)); %Plummer & Busenberg 1982 (10.1016/0016-7037(82)90056-4)
    T.Kh_CH4=exp(-68.8862+101.4956.*100./(T.Temp_Kelvin)+28.7314.*log((T.Temp_Kelvin)./100))./((T.Temp_Kelvin).*gas_constant./100); %Wiesenburg & Guinasso 1979 (10.1021/je60083a006) (taken from Wanninkhof 2014)
    T.Kh_N2O=exp(-62.7062+97.3066.*100./(T.Temp_Kelvin)+24.1406.*log((T.Temp_Kelvin)./100))./((T.Temp_Kelvin).*gas_constant./100); %Weiss & Price 1980 (10.1016/0304-4203(80)90024-9) (taken from Wanninkhof 2014)
    for m=1:size(T(:,1),1)
        T.CO2_atm(m)=CO2_sat.*T.Kh_CO2(m);
        T.CH4_atm(m)=CH4_sat.*T.Kh_CH4(m);
        T.N2O_atm(m)=N2O_sat.*T.Kh_N2O(m);
    end

end



%% 1. Find matching latitude and longitude and large river sites that are not an exact lat/long match, but less than 5km apart
% For lakes/reservoirs, simply merge exact matches for now (merging by lake is done through another function)

if isequal(choice,'streams/rivers')
    %find sites with inconsistent stream order values and replace
    [~,~,idx_site]=unique([T.Latitude,T.Longitude],'rows');
    n_fixed=0;
    for site_id=1:max(idx_site)
        rows=find(idx_site==site_id);
        if numel(unique(T.streamorder_consolidated(rows)))>1
            strahler_vals=unique(T.StrahlerOrder(rows));
            strahler_vals(isnan(strahler_vals))=[]; %remove NaNs
            if numel(strahler_vals)==1
                T.streamorder_consolidated(rows)=strahler_vals;
                n_fixed=n_fixed+numel(rows);
            end
        end
    end

    %nearby location merge for large rivers
    coords=[T.Latitude,T.Longitude];
    large_river_mask=T.streamorder_consolidated>SO_large;
    coords_large_rivers=coords(large_river_mask,:);
    n=size(coords_large_rivers,1);
    %compute pairwise Haversine distance matrix
    D=zeros(n,n);
    for i=1:n
        for j=i+1:n
            D(i,j)=lldistkm(coords_large_rivers(i,:),coords_large_rivers(j,:));
            D(j,i)=D(i,j);
        end
    end
    %merge with dbscan
    [labels_large,~]=dbscan(D,threshold_km,minpts,'Distance','precomputed');
    labels_large(labels_large>0)=labels_large(labels_large>0)+100000; %add 10000 to make sure we don't get any overlap
    labels_large(labels_large==-1)=0;
    idx_combined=zeros(height(T),1);
    idx_combined(large_river_mask)=labels_large;
    unlabeled_mask=idx_combined==0;
    coords_unlabeled=coords(unlabeled_mask,:);

    %exact location merge for smaller rivers
    [~,~,exact_labels]=unique(coords_unlabeled,'rows');
    idx_combined(unlabeled_mask)=exact_labels+200000; %add 20000 to make sure we don't get any overlap

    %exact location merge for all dataset to indirectly determine the number of nearby merges
    Loc=[T.Latitude,T.Longitude];
    [~,~,idx_exact]=unique(Loc,'rows','stable');

    %display results
    fprintf('Initial dataset size: %d\n',height(T));
    fprintf('Unique site count after exact merging: %d\n',numel(unique(idx_exact)));
    fprintf('Unique site count after exact and nearby merging: %d\n',numel(unique(idx_combined)));
    fprintf('Total duplicates removed: %d\n',height(T)-numel(unique(idx_combined)));
    pause(0.1)


    %plot map
    figure;
    geoscatter(T.Latitude,T.Longitude,50,[.5 .5 .5],'filled');
    hold on;
    exact_counts=accumarray(idx_exact,1);
    exact_repeats=ismember(idx_exact,find(exact_counts>1));
    nearby_merge_mask=idx_combined>=100001 & idx_combined<200000;
    geoscatter(T.Latitude(nearby_merge_mask),T.Longitude(nearby_merge_mask),51,[1 .6 .6],'filled')%,'MarkerFaceAlpha',0.5);
    geoscatter(T.Latitude(exact_repeats),T.Longitude(exact_repeats),52,[.5 .5 .5],'filled');
    legend({'all sites','merged sites'},'Location','best','FontSize',9);
    geobasemap('topographic');

    % Simple merge by lake index for lakes/reservoirs
elseif isequal(choice,'lakes/reservoirs')
    [~,~,idx_combined]=unique(T.Lake_idx);
end



%% 2. Other pre-calculations
[~,unique_idx]=unique(idx_combined,'stable'); %get indices
T.SiteID=idx_combined;

%compute dates as datenum
T.date=datenum(T.year,T.month,T.day);

%get numerical values for direct/indirect CO2 measurements, and apply one value per site as follows:
% if all measurements are direct, enter 1 (direct); if at least one measurement is indirect, enter 0 (indirect)
T.CO2_method_num=NaN(size(T,1),1); %preallocate
T.CO2_method_num(T.CO2_method_cat=='direct')=1;
T.CO2_method_num(T.CO2_method_cat=='indirect')=0;
T=movevars(T,"CO2_method_num",'After',"CO2_method_cat");
for ll=1:max(unique_idx)
    E=T.CO2_method_num(idx_combined==ll);
    if sum(E)<size(E,1)
        T.CO2_method_num(idx_combined==ll)=0; %here we are conservative; if not all measurements are direct, we consider them as indirect
    elseif sum(E)==size(E,1)
        T.CO2_method_num(idx_combined==ll)=1;
    end
end


%% 3. Calculate means/medians of all measurements for single sites and for a range of variables

if isequal(choice,'streams/rivers')
    R=[T.Latitude T.Longitude T.date T.Discharge_converted T.CO2_converted T.F_CO2_converted T.CH4_converted...
        T.F_CH4_converted T.F_CH4_ebb_converted T.N2O_converted T.F_N2O_converted...
        T.k600_converted T.Climate T.streamorder_consolidated T.slope_consolidated T.RiverWidth T.WatershedArea T.Temperature...
        T.CO2_method_num T.humanFootprint_BasAt T.soilOrganicCarbon_BasAt T.airTemperature_BasAt...
        T.rainfall_BasAt T.runoff_BasAt T.wetlandExtent_BasAt T.forestExtent_BasAt...
        T.aridity_BasAt T.groundwaterDepth_BasAt T.netPrimaryProductivity_Modis...
        T.peatlandCover_Modis T.elevation_H90 T.catchmentArea_H90...
        T.catchmentArea_GRADES T.discharge T.velocity T.Seasonality T.DO_converted T.pH...
        T.Conductivity T.DOC_converted T.DIC_converted];
elseif isequal(choice,'lakes/reservoirs')
    R=[T.Latitude T.Longitude T.date T.CO2_converted T.F_CO2_converted T.CH4_converted T.F_CH4_converted...
        T.F_CH4_ebb_converted T.N2O_converted T.F_N2O_converted T.CO2_atm T.CH4_atm T.N2O_atm...
        T.k600_converted T.k600_consolidated T.kco2_W92 T.kch4_W92 T.kn2o_W92 T.kco2_CC98 T.kch4_CC98 T.kn2o_CC98...
        T.kco2_VP13 T.kch4_VP13 T.kn2o_VP13 T.CO2_method_num T.Climate...
        T.Temperature T.area_consolidated T.MeanDepth T.MaxDepth T.Volume,T.DO_converted,T.pH,T.Conductivity];
end

Av=NaN(size(unique_idx,1),size(R,2)); %preallocate
As=NaN(size(unique_idx,1),size(R,2));
[~,~,idx_new]=unique(idx_combined);

meanmed=questdlg('Use mean or median?','averaging','mean','median','median');
switch meanmed
    case 'mean'
        for d=1:size(R,2)
            Av(:,d)=accumarray(idx_new,R(:,d),[size(unique_idx,1) 1],@(x)mean(x,'omitnan')); %calculate means while removing NaNs
            As(:,d)=accumarray(idx_new,R(:,d),[size(unique_idx,1) 1],@(x)sum(~isnan(x))); %calculate number of elements for each mean calculation
        end
    case 'median'
        for d=1:size(R,2)
            Av(:,d)=accumarray(idx_new,R(:,d),[size(unique_idx,1) 1],@(x)median(x,'omitnan')); %calculate medians while removing NaNs
            As(:,d)=accumarray(idx_new,R(:,d),[size(unique_idx,1) 1],@(x)sum(~isnan(x))); %calculate number of elements for each median calculation
        end
end

% extract the measurement/reported frequencies for each site, and pick the highest frequency when sites have different frequencies across studies
freq_order=["HF","subdaily","daily","weekly","fortnightly","monthly","seasonal","annual","interannual","single"];
freq_rank_map=containers.Map(freq_order,1:length(freq_order));
freq_meas=string(T.measurement_frequency);
freq_rep=string(T.reported_frequency);
all_freq_meas=splitapply(@(x) strjoin(unique(x),'; '),freq_meas,idx_new);
all_freq_rep=splitapply(@(x) strjoin(unique(x),'; '),freq_rep,idx_new);
highest_freq_meas=splitapply(@(x) getHighestFreq(x,freq_rank_map),freq_meas,idx_new);
highest_freq_rep=splitapply(@(x) getHighestFreq(x,freq_rank_map),freq_rep,idx_new);
% extract flux method for each site, and pick 'measured' if at least one 'measured' flux.
F_method_site=splitapply(@(x) getFluxMethod(x),T.F_method_cat,idx_new);

%extract other variables for merged sites, pick most common (KG) and concatenate (Ref and SiteName)
KG_site=splitapply(@(x) getMostCommon(x),string(T.KG),idx_new);
Ref_site=splitapply(@(x) strjoin(unique(x),'; '),string(T.Ref),idx_new);
SiteName_site=splitapply(@(x) strjoin(unique(x),'; '),string(T.SiteName),idx_new);

if isequal(choice,'lakes/reservoirs')
    %extract system type for lakes (most common category) and lake index (most common but should be always the same)
    LakeType_site=splitapply(@(x) getMostCommon(x),T.LakeType,idx_new);
    LakeIdx_site=splitapply(@(x) getMostCommon(x),T.Lake_idx,idx_new);
end

%% 4. Create a smaller table with averaged variables
if isequal(choice,'streams/rivers')
    vars={'Lat','Long','Mid_date','Discharge_converted','CO2_converted','F_CO2_converted',...
        'CH4_converted','F_CH4_diff_converted','F_CH4_ebb_converted','N2O_converted',...
        'F_N2O_converted','k600_converted','Climate','streamorder_consolidated','slope_consolidated','Width','Catchment_area','Temperature',...
        'CO2_direct_indirect','humanFootprint_BasAt','soilOrganicCarbon_BasAt',...
        'airTemperature_BasAt','rainfall_BasAt','runoff_BasAt','wetlandExtent_BasAt',...
        'forestExtent_BasAt','aridity_BasAt','groundwaterDepth_BasAt',...
        'netPrimaryProductivity_Modis','peatlandCover_Modis',...
        'elevation_H90','catchmentArea_H90',...
        'catchmentArea_GRADES','discharge','velocity','seasonality','DO_converted','pH',...
        'Conductivity','DOC_converted','DIC_converted'};
    nbs={'Nb_CO2','Nb_FCO2','Nb_CH4','Nb_FCH4diff','Nb_FCH4ebb','Nb_N2O','Nb_FN2O'};
    nb_indices=[5,6,7,8,9,10,11];
    Av_table=array2table(Av,'VariableNames',vars);
    As_table=array2table(As(:,nb_indices),'VariableNames',nbs);
    Ts=[Av_table,As_table];
    Ts.measurement_frequency_combined=all_freq_meas;
    Ts.highest_measurement_frequency=highest_freq_meas;
    Ts.reported_frequency_combined=all_freq_rep;
    Ts.highest_reported_frequency=highest_freq_rep;
    Ts.KG=KG_site;
    Ts.Ref=Ref_site;
    Ts.SiteName=SiteName_site;
    Ts.F_method_cat=F_method_site;
    Ts=movevars(Ts,"Ref",'Before',"Lat");
    Ts=movevars(Ts,"SiteName",'After',"Ref");
    Ts=movevars(Ts,"KG",'After',"Mid_date");
    Ts=movevars(Ts,"F_method_cat",'After',"KG");

elseif isequal(choice,'lakes/reservoirs')       
          vars={'Lat','Long','Mid_date','CO2_converted','F_CO2_converted',...
        'CH4_converted','F_CH4_diff_converted','F_CH4_ebb_converted',...
        'N2O_converted','F_N2O_converted','CO2_atm','CH4_atm','N2O_atm','k600_converted',...
        'k600_consolidated','kCO2_W92','kCH4_W92','kN2O_W92','kCO2_CC98','kCH4_CC98','kN2O_CC98'...
        'kCO2_VP13','kCH4_VP13','kN2O_VP13','CO2_method_num','Climate','Temperature',...
        'Area','MeanDepth','MaxDepth','Volume','DO_converted','pH','Conductivity'};
    nbs={'Nb_CO2','Nb_FCO2','Nb_CH4','Nb_FCH4diff','Nb_FCH4ebb','Nb_N2O','Nb_FN2O'};
    nb_indices=[4,5,6,7,8,9,10];
    Av_table=array2table(Av,'VariableNames',vars);
    As_table=array2table(As(:,nb_indices),'VariableNames',nbs);
    Ts=[Av_table,As_table];
    Ts.measurement_frequency_combined=all_freq_meas;
    Ts.highest_measurement_frequency=highest_freq_meas;
    Ts.reported_frequency_combined=all_freq_rep;
    Ts.highest_reported_frequency=highest_freq_rep;
    Ts.KG=KG_site;
    Ts.Ref=Ref_site;
    Ts.SiteName=SiteName_site;
    Ts.F_method_cat=F_method_site;
    Ts.LakeType=LakeType_site;
    Ts.Lake_idx=LakeIdx_site;
    Ts=movevars(Ts,"Ref",'Before',"Lat");
    Ts=movevars(Ts,"SiteName",'After',"Ref");
    Ts=movevars(Ts,"KG",'After',"Mid_date");
    Ts=movevars(Ts,"F_method_cat",'After',"KG");
    Ts=movevars(Ts,"LakeType",'After',"F_method_cat");
    Ts=movevars(Ts,"Lake_idx",'After',"SiteName");
end

% trim off dataset to keep only (sub)tropics
Ts=Ts(Ts.Lat<34.00001 & Ts.Lat>-34.00001,:);

%output text file
if isequal(choice,'lakes/reservoirs')
    txtout='Clean_data_lakes_reservoirs_reduced.txt';
    txtfile=fopen(txtout,'w');
    fprintf(txtfile,'ref\tsite\tdate\tlatitude\tlongitude\tlake_type\tclimate_class\tarea_km2\tmean_depth_m\n');
    for u=1:size(Ts,1)
        fprintf(txtfile,'%s\t%s\t%f\t%f\t%f\t%s\t%f\t%f\t%f\n',Ts.Ref{u},Ts.SiteName{u},Ts.Mid_date(u),Ts.Lat(u),Ts.Long(u),...
            Ts.LakeType{u},Ts.Climate(u),Ts.Area(u),Ts.MeanDepth(u));
    end
    fclose(txtfile);
end


disp(' ')
disp('calculation of means/medians per site done')
disp(' ')

end


