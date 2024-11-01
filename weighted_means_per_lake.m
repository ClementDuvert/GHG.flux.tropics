%% Function to calculate weighted mean values for each lake/reservoir with several measurements
%
function[T]=weighted_means_per_lake()
%
% Inputs
% Requires the following files in the 'input datasets' folder:
% 'Clean_data_lakes_input.txt' (clean file with added HydroLAKES data)
% 'Lake_properties_manual.csv' (surface area estimates for systems not in HydroLAKES)
% 'Lake_wind.csv'              (WorldClim2 wind data for each system)
%
% Outputs
% T is a table with one row per lake/reservoir and the respective weighted mean values



%% 1. Input datasets + initial calculations

disp(' '); disp('starting calculation of weighted means per system...'); disp(' ')

%constants
CO2_sat=415; %umol/mol or ppmv
CH4_sat=1.87; %umol/mol or ppmv
N2O_sat=0.332; %umol/mol or ppmv  in 2020
gas_constant=8.3144626181; %m3 Pa mol-1 K-1

% import dataset with HydroLAKES parameters included
SA=readtable('Clean_data_lakes_input.txt');

% import dataset for lakes that are not in HydroLAKES
SB=readtable('Lake_properties_manual.csv');

%import wind speed data and incorporate into SA table
SC=readtable('Lake_wind.csv');
SA.wind=NaN(size(SA,1),1);
for u=1:size(SA,1)
    matchIdx=SC.lake_ID==SA.group_idx_manual(u);
    if any(matchIdx)
        SA.wind(u)=mean(SC.wind_worldclim(matchIdx));
        clear matchIdx
    end
end

%gap fill temperature data if no data
SA.Temp_gapfill=SA.temperature;
a=repelem(20,size(SA,1)); %assign temperature to data with no reported temperature
b=isnan(SA.Temp_gapfill); SA.Temp_gapfill(b)=a(b);
% calculate Schmidt numbers
SA.Sc_CO2=1923.6-125.06.*SA.Temp_gapfill+4.3773.*SA.Temp_gapfill.^2-0.085681.*SA.Temp_gapfill.^3+0.00070284.*SA.Temp_gapfill.^4; %Jahne et al (1987) / Wanninkhof (2014)
SA.Sc_CH4=1909.4-120.78.*SA.Temp_gapfill+4.1555.*SA.Temp_gapfill.^2-0.080578.*SA.Temp_gapfill.^3+0.00065777.*SA.Temp_gapfill.^4; %Jahne et al (1987) / Wanninkhof (2014)
SA.Sc_N2O=2141.2-152.56.*SA.Temp_gapfill+5.8963.*SA.Temp_gapfill.^2-0.12411.*SA.Temp_gapfill.^3+0.0010655.*SA.Temp_gapfill.^4; %Jahne et al (1987) / Wanninkhof (2014)

% calculate k600 based on Wanninkhof 1992; Cole and Caraco 1998; Vachon and Prairie 2013
%this requires adding lake areas first
SA.size=SA.hylak_lake_area;
for s=1:size(SA,1)
    if isnan(SA.size(s))
        SA.size(s)=SA.surface_area_km2(s);
        if isnan(SA.surface_area_km2(s))
            SA.size(s)=SB.lake_area_km2(SB.group_idx_manual==SA.group_idx_manual(s));
        end
    end
end
%then we compute k600
for a=1:size(SA,1)
    w1=0.45*SA.wind(a)^1.64; %in cm/h
    SA.k600_modelled_W92(a)=w1/100*24; %in m/d
    w2=2.07+0.215*SA.wind(a)^1.7; %in cm/h
    SA.k600_modelled_CC98(a)=w2/100*24; %in m/d
    w3=2.51+(1.48*SA.wind(a))+(0.39*SA.wind(a)*log10(SA.size(a))); %in cm/h
    SA.k600_modelled_VP13(a)=w3/100*24; %in m/d
    clear w1 w2 w3
end

% calculate kCO2, kCH4 and kN2O
gases={'CO2','CH4','N2O'};
models={'W92','CC98','VP13'};
for i=1:length(models)
    for j=1:length(gases)
        gas=gases{j};
        model=models{i};
        SA.(['k' lower(gas) '_' model])=SA.(['k600_modelled_' model]).*(SA.(['Sc_' gas])/600).^-0.5;
    end
end

% calculate Henry constants and atmospheric concentrations
SA.Temp_Kelvin=273.15+SA.Temp_gapfill;
SA.Kh_CO2=10.^-(-(108.3865+0.01985076.*(SA.Temp_Kelvin)-6919.53./(SA.Temp_Kelvin)-40.4515.*log10(SA.Temp_Kelvin)+669365./(SA.Temp_Kelvin).^2)); %Plummer & Busenberg 1982 (10.1016/0016-7037(82)90056-4)
SA.Kh_CH4=exp(-68.8862+101.4956.*100./(SA.Temp_Kelvin)+28.7314.*log((SA.Temp_Kelvin)./100))./((SA.Temp_Kelvin).*gas_constant./100); %Wiesenburg & Guinasso 1979 (10.1021/je60083a006) (taken from Wanninkhof 2014)
SA.Kh_N2O=exp(-62.7062+97.3066.*100./(SA.Temp_Kelvin)+24.1406.*log((SA.Temp_Kelvin)./100))./((SA.Temp_Kelvin).*gas_constant./100); %Weiss & Price 1980 (10.1016/0304-4203(80)90024-9) (taken from Wanninkhof 2014)
for m=1:size(SA(:,1),1)
    SA.CO2_atm(m)=CO2_sat.*SA.Kh_CO2(m);
    SA.CH4_atm(m)=CH4_sat.*SA.Kh_CH4(m);
    SA.N2O_atm(m)=N2O_sat.*SA.Kh_N2O(m);
end


%% 2. Calculate weighted mean concentrations and fluxes for each lake/reservoir
s=unique(SA.group_idx_manual); %find individual lake IDs
e=1;
n_lakes=numel(s); %preallocate
T=NaN(n_lakes,35); 
lake_refs=cell(n_lakes,1);

variables={'CO2_umolL_1','CH4_umolL_1','N2O_umolL_1','FCO2_mmolm_2d_1','FCH4diff_mmolm_2d_1','FCH4ebb_mmolm_2d_1','FN2O_mmolm_2d_1','CO2_atm','CH4_atm','N2O_atm'};
k_variables={'kco2_W92','kch4_W92','kn2o_W92','kco2_CC98','kch4_CC98','kn2o_CC98','kco2_VP13','kch4_VP13','kn2o_VP13'}; 
for x=s(1):s(end)
    Cut=SA(SA.group_idx_manual==x,:); %trim dataset for lake of interest
    if ~isempty(Cut)
        weighted_vars=NaN(1,length(variables)); %preallocate
        nb_meas_tot=NaN(1,length(variables));

        % calculate weighted means for each variable
        for v=1:length(variables)
            weights=Cut.(['nb_meas_' extractBefore(variables{v},'_')])(~isnan(Cut.(variables{v})));
            array=Cut.(variables{v})(~isnan(Cut.(variables{v})));
            if ~isempty(weights)
                weighted_vars(v)=sum(weights.*array)/sum(weights);
                nb_meas_tot(v)=sum(weights);
            end
        end

        weighted_k_vars=NaN(1,length(k_variables)); %preallocate
        nb_meas_k_tot=NaN(1,length(k_variables)); 
        for v=1:length(k_variables)
            k_var_name=k_variables{v};
            if contains(k_var_name,'co2','IgnoreCase',true)
                gas='CO2';
            elseif contains(k_var_name,'ch4','IgnoreCase',true)
                gas='CH4';
            elseif contains(k_var_name,'n2o','IgnoreCase',true)
                gas='N2O';
            end
            weights=Cut.(['nb_meas_' gas])(~isnan(Cut.(k_var_name)));
            array=Cut.(k_var_name)(~isnan(Cut.(k_var_name)));
            if ~isempty(weights)
                weighted_k_vars(v)=sum(weights.* array)/sum(weights);
                nb_meas_k_tot(v)=sum(weights);
            end
        end
      
        % extract lake information
        lake_ID=Cut.group_idx_manual(1);
        lake_climate=Cut.climate_class(1);
        lake_latitude=Cut.latitude(1);
        lake_longitude=Cut.longitude(1);
        lake_area=Cut.size(1);
        lake_site=Cut.site(1);
        lake_type=Cut.lake_type(1);
        lake_k600=Cut.k600_md_1(1);
        lake_depth_hylak=Cut.hylak_depth_avg(1);
        lake_depth_ref=Cut.mean_depth_m(1);
        lake_ref=strjoin(arrayfun(@num2str,unique(Cut.ref),'UniformOutput',false),';'); %we include all references where data were extracted from
        
        %compute CO2 method (direct/indirect)
        CO2_method=NaN;
        s1=mean(Cut.CO2_method(~isnan(Cut.CO2_method)));
        if s1>0.5
            CO2_method=1;
        elseif s1<=0.5
            CO2_method=0;
        end
        clear s1

        % concatenate and add sites, types, refs        
        T(e,:)=[lake_ID,lake_climate,lake_depth_hylak,lake_depth_ref,lake_latitude,lake_longitude,...
            lake_area,weighted_vars,weighted_k_vars,nb_meas_tot(:,1:7),lake_k600,CO2_method];
        Sites{e}=lake_site;
        Types{e}=lake_type;
        lake_refs{e}=lake_ref;
        e=e+1;
    end
end

% convert to table
T=array2table(T,'VariableNames',{'lake_ID','climate','depth_hylak','depth_ref',...
    'latitude','longitude','lake_area','weighted_pCO2','weighted_pCH4','weighted_pN2O','weighted_FCO2',...
    'weighted_FCH4_diff','weighted_FCH4_ebb','weighted_FN2O','weighted_CO2atm','weighted_CH4atm',...
    'weighted_N2Oatm','weighted_kCO2_W92','weighted_kCH4_W92','weighted_kN2O_W92',...
    'weighted_kCO2_CC98','weighted_kCH4_CC98','weighted_kN2O_CC98','weighted_kCO2_VP13',...
    'weighted_kCH4_VP13','weighted_kN2O_VP13','nb_meas_tot_CO2','nb_meas_tot_CH4',...
    'nb_meas_tot_N2O','nb_meas_tot_FCO2','nb_meas_tot_FCH4diff','nb_meas_tot_FCH4ebb',...
    'nb_meas_tot_FN2O','k600','CO2_method'});

T.ref=lake_refs(1:e-1);
for x=1:size(T,1)
    T.site(x)=cellstr(Sites{x});
    T.type(x)=cellstr(Types{x});
end


%% 3. Export data
% option to add modelled fluxes based on k and concentrations, or export DB with measurements only
choice=questdlg('Model fluxes?','Fluxes','only measurements','include modelled fluxes','include modelled fluxes');
switch choice

    case 'only measurements'
        %export dataset
        cd('/Users/ClementDuvert/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/Documents/CDU/RESEARCH/PROJECTS/12. DECRA/Global Database/MATLAB/output_datasets')
        txtout='Clean_data_lakes_reservoirs_per_lake_ALL_onlymeasuredF.txt';
        txtfile=fopen(txtout,'w');
        fprintf(txtfile,'ref\tsite\tlatitude\tlongitude\tclimate_class\tlake_ID\tlake_area_km2\tlake_depth_hylak_m\tlake_depth_ref_m\tlake_type\tCO2_umolL-1\tnb_meas_CO2\tCO2_method\tFCO2_mmolm-2d-1\tnb_meas_FCO2\tCH4_umolL-1\tnb_meas_CH4\tFCH4diff_mmolm-2d-1\tnb_meas_FCH4diff\tFCH4ebb_mmolm-2d-1\tnb_meas_FCH4ebb\tN2O_umolL-1\tnb_meas_N2O\tFN2O_mmolm-2d-1\tnb_meas_FN2O\n');
        for u=1:size(T,1)
            fprintf(txtfile,'%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',T.ref{u},T.site{u},T.latitude(u),T.longitude(u),T.climate(u),T.lake_ID(u),T.lake_area(u),T.depth_hylak(u),T.depth_ref(u),T.type{u},T.weighted_pCO2(u),T.nb_meas_tot_CO2(u),T.CO2_method(u),T.weighted_FCO2(u),T.nb_meas_tot_FCO2(u),T.weighted_pCH4(u),T.nb_meas_tot_CH4(u),T.weighted_FCH4_diff(u),T.nb_meas_tot_FCH4diff(u),T.weighted_FCH4_ebb(u),T.nb_meas_tot_FCH4ebb(u),T.weighted_pN2O(u),T.nb_meas_tot_N2O(u),T.weighted_FN2O(u),T.nb_meas_tot_FN2O(u));
        end
        fclose(txtfile);

    case 'include modelled fluxes'
        % pre-assign columns for each model
        T.weighted_FCO2_W92=NaN(size(T,1),1); T.weighted_FCO2_CC98=NaN(size(T,1),1); T.weighted_FCO2_VP13=NaN(size(T,1),1);
        T.weighted_FCH4_diff_W92=NaN(size(T,1),1); T.weighted_FCH4_diff_CC98=NaN(size(T,1),1); T.weighted_FCH4_diff_VP13=NaN(size(T,1),1);
        T.weighted_FN2O_W92=NaN(size(T,1),1); T.weighted_FN2O_CC98=NaN(size(T,1),1); T.weighted_FN2O_VP13=NaN(size(T,1),1);

        % calculate and add mean fluxes for lakes where only concentrations are available
        for lk=1:size(T,1)
            if isnan(T.weighted_FCO2(lk)) && ~isnan(T.weighted_pCO2(lk))
                T.weighted_FCO2_W92(lk)=(T.weighted_pCO2(lk)-T.weighted_CO2atm(lk))*T.weighted_kCO2_W92(lk); %uM * m/d = mmol/m2/d
                T.weighted_FCO2_CC98(lk)=(T.weighted_pCO2(lk)-T.weighted_CO2atm(lk))*T.weighted_kCO2_CC98(lk); %uM * m/d = mmol/m2/d
                T.weighted_FCO2_VP13(lk)=(T.weighted_pCO2(lk)-T.weighted_CO2atm(lk))*T.weighted_kCO2_VP13(lk); %uM * m/d = mmol/m2/d
                T.weighted_FCO2(lk)=mean([T.weighted_FCO2_CC98(lk),T.weighted_FCO2_VP13(lk)]); %we take the mean of the best 2 models and add to flux arrays
            end
            if isnan(T.weighted_FCH4_diff(lk)) && ~isnan(T.weighted_pCH4(lk))
                T.weighted_FCH4_diff_W92(lk)=(T.weighted_pCH4(lk)-T.weighted_CH4atm(lk))*T.weighted_kCH4_W92(lk);
                T.weighted_FCH4_diff_CC98(lk)=(T.weighted_pCH4(lk)-T.weighted_CH4atm(lk))*T.weighted_kCH4_CC98(lk);
                T.weighted_FCH4_diff_VP13(lk)=(T.weighted_pCH4(lk)-T.weighted_CH4atm(lk))*T.weighted_kCH4_VP13(lk);
                T.weighted_FCH4_diff(lk)=mean([T.weighted_FCH4_diff_CC98(lk),T.weighted_FCH4_diff_VP13(lk)]); %we take the mean of the best 2 models and add to flux arrays
            end
            if isnan(T.weighted_FN2O(lk)) && ~isnan(T.weighted_pN2O(lk))
                T.weighted_FN2O_W92(lk)=(T.weighted_pN2O(lk)-T.weighted_N2Oatm(lk))*T.weighted_kN2O_W92(lk);
                T.weighted_FN2O_CC98(lk)=(T.weighted_pN2O(lk)-T.weighted_N2Oatm(lk))*T.weighted_kN2O_CC98(lk);
                T.weighted_FN2O_VP13(lk)=(T.weighted_pN2O(lk)-T.weighted_N2Oatm(lk))*T.weighted_kN2O_VP13(lk);
                T.weighted_FN2O(lk)=mean([T.weighted_FN2O_CC98(lk),T.weighted_FN2O_VP13(lk)]); %we take the mean of the best 2 models and add to flux arrays
            end
        end

        %export dataset
        %cd('/Users/ClementDuvert/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/Documents/CDU/RESEARCH/PROJECTS/12. DECRA/Global Database/MATLAB/output_datasets')
        txtout='Clean_data_lakes_reservoirs_per_lake_ALL_withmodelledF.txt';
        txtfile=fopen(txtout,'w');
        fprintf(txtfile,'ref\tsite\tlatitude\tlongitude\tclimate_class\tlake_ID\tlake_area_km2\tlake_depth_hylak_m\tlake_depth_ref_m\tlake_type\tCO2_umolL-1\tnb_meas_CO2\tCO2_method\tFCO2_mmolm-2d-1\tFCO2_mmolm-2d-1_W92\tFCO2_mmolm-2d-1_CC98\tFCO2_mmolm-2d-1_VP13\tnb_meas_FCO2\tCH4_umolL-1\tnb_meas_CH4\tFCH4diff_mmolm-2d-1\tFCH4diff_mmolm-2d-1_W92\tFCH4diff_mmolm-2d-1_CC98\tFCH4diff_mmolm-2d-1_VP13\tnb_meas_FCH4diff\tFCH4ebb_mmolm-2d-1\tnb_meas_FCH4ebb\tN2O_umolL-1\tnb_meas_N2O\tFN2O_mmolm-2d-1\tFN2O_mmolm-2d-1_W92\tFN2O_mmolm-2d-1_CC98\tFN2O_mmolm-2d-1_VP13\tnb_meas_FN2O\n');
        for u=1:size(T,1)
            fprintf(txtfile,'%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',T.ref{u},T.site{u},T.latitude(u),T.longitude(u),T.climate(u),T.lake_ID(u),T.lake_area(u),T.depth_hylak(u),T.depth_ref(u),T.type{u},T.weighted_pCO2(u),T.nb_meas_tot_CO2(u),T.CO2_method(u),T.weighted_FCO2(u),T.weighted_FCO2_W92(u),T.weighted_FCO2_CC98(u),T.weighted_FCO2_VP13(u),T.nb_meas_tot_FCO2(u),T.weighted_pCH4(u),T.nb_meas_tot_CH4(u),T.weighted_FCH4_diff(u),T.weighted_FCH4_diff_W92(u),T.weighted_FCH4_diff_CC98(u),T.weighted_FCH4_diff_VP13(u),T.nb_meas_tot_FCH4diff(u),T.weighted_FCH4_ebb(u),T.nb_meas_tot_FCH4ebb(u),T.weighted_pN2O(u),T.nb_meas_tot_N2O(u),T.weighted_FN2O(u),T.weighted_FN2O_W92(u),T.weighted_FN2O_CC98(u),T.weighted_FN2O_VP13(u),T.nb_meas_tot_FN2O(u));
        end
        fclose(txtfile);

        disp(' ')
        disp('calculation of weighted means per lake done')
        disp(' ')

end
