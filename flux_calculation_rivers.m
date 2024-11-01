%% Function to upscale GHG fluxes per stream order and climate zone using a bootstrapping approach
%
function[T,arealFlux_CO2,scaledFlux_CO2,scaledFlux_eq_CO2,arealFlux_CH4,...
    scaledFlux_CH4,scaledFlux_eq_CH4,arealFlux_N2O,scaledFlux_N2O,...
    scaledFlux_eq_N2O]=flux_calculation_rivers(T)
%
% This function contains the following sections:
% 1. Initialise constants and import datasets
% 2. Preliminary calculations (Sc, k, atmospheric gas)
% 3. Estimate fluxes based on concentration and k
% 4. Bootstrapping to obtain median (Q1-Q3) fluxes per stream order and climate zone
% 5. Upscaling using river surface areas
%
% Inputs
% T is a reduced table that has been processed using 'unit_conversion', 'climate_classification' and 'medians_per_site'
% In addition, the following files are needed in the input dataset folder:
% 'Hydro_properties_BasinATLAS.csv'
% 'Hydro_properties_npp_peat.csv'
% 'Hydro_properties_H90.csv'
% 'Hydro_properties_GRADES.csv'
% 'river_surface_areas.csv'
%
% Outputs
% all outputs are organised as follows. 3 rows: Q1, med, Q3. 9 columns: stream order 1 to 9+. 5 depths: climate zones.
% arealFlux_gas is a 3x9x5 array with areal fluxes in mmol/m2/d. 
% scaledFlux_gas is is a 3x9x5 array with upscaled fluxes in Tg C/yr or Tg N/yr.
% scaledFlux_eq_gas is a 3x9x5 array with upscaled fluxes in Tg CO2-eq/yr.


%% 1. Initialise constants and import datasets

disp(' '); disp('starting river flux calculation...'); disp(' ')

%constants
boot_nb=1000; %nb of bootstrap samples
gravit=9.80665; %m/s2
molar_mass_C=12.011; %g mol-1
molar_mass_N=14.01; %g mol-1
molar_mass_CO2=44.01; %g mol-1
molar_mass_CH4=16.05; %g mol-1
molar_mass_N2O=44.02; %g mol-1
CO2_sat=415; %umol/mol or ppmv
CH4_sat=1.87; %umol/mol or ppmv
N2O_sat=0.332; %umol/mol or ppmv  in 2020
gas_constant=8.3144626181; %m3 Pa mol-1 K-1
GWP_N2O=298;
GWP_CH4=34;
percent_error_area=30; %percent error on surface area estimates
climate_class={'humid tropics','wet-dry tropics','(semi)arid (sub)tropics','humid subtropics','highland (sub)tropics'};

% create row with site indices to match indices in BasinATLAS file
T.Site_ID=(1:1:size(T,1))';
T=movevars(T,"Site_ID",'Before',"Ref");

% import and add landscape properties (BasinATLAS)
BasAt=readtable('Hydro_properties_BasinATLAS.csv');
BasAt=sortrows(BasAt,'site_numeric'); %sort rows according to site indices
T.humanFootprint_BasAt=BasAt.human_footprint_up_09;
T.soilOrganicCarbon_BasAt=BasAt.soil_org_carbon_up_tonneshectare;
T.airTemperature_BasAt=BasAt.air_temp_up_avg_celsius/10;
T.rainfall_BasAt=BasAt.precip_up_mm;
T.runoff_BasAt=BasAt.runoff_mm;
T.wetlandExtent_BasAt=BasAt.wetland_grouped2_up_per;
T.forestExtent_BasAt=BasAt.forest_cover_up_per;
T.aridity_BasAt=BasAt.aridity_index_up;
T.groundwaterDepth_BasAt=BasAt.gw_table_cm;

% import and add landscape properties (Modis)
Modis=readtable('Hydro_properties_npp_peat.csv');
Modis=sortrows(Modis,'site_numeric'); %sort rows according to site indices
T.netPrimaryProductivity_Modis=Modis.NPP_yr;
T.peatlandCover_Modis=Modis.peatland_cover;

% import and add hydro properties (Hydrography90m)
H90=readtable('Hydro_properties_H90.csv');
H90=sortrows(H90,'site_numeric'); %sort rows according to site indices
T.slope_H90=H90.slope_channel;
T.elevation_H90=H90.elevation_m;
T.streamOrder_H90=H90.order_snapped_coords;
T.shifted_streamOrder_H90=T.streamOrder_H90-1;
T.shifted_streamOrder_H90(~isnan(T.StrahlerOrder))=T.StrahlerOrder(~isnan(T.StrahlerOrder))-1; %correct stream orders with actual order if known
% 'T.shifted_streamOrder_H90' is the variable we use from now on as the correct SOs
T.catchmentArea_H90=H90.catchment_area_snapped; %km2

% import other properties from GRADES
GRADES=readtable('Hydro_properties_GRADES.csv');
GRADES=sortrows(GRADES,'site_numeric'); %sort rows according to site indices
T.slope_GRADES=GRADES.slope_grades;
T.slope_GRADES(T.slope_GRADES==0)=5E-5; %assign very low values to sites with slope 0
T.slope_GRADES(~isnan(T.Slope))=T.Slope(~isnan(T.Slope)); %correct slopes with actual value if known
T.catchmentArea_GRADES=GRADES.upstream_area_grades/1E4;
T.streamOrder_GRADES=GRADES.order_grades;

%import river surface areas for each SO and climate class
SurfAr=readtable('river_surface_areas.csv');

% calculate discharge and flow velocity based on BasinATLAS (runoff) and H90 (area)
T.discharge=(T.runoff_BasAt/1E3).*(T.catchmentArea_H90*1E6)/(3600*24*365); %discharge in m3/s
T.velocity=exp(0.12*log(T.discharge)-1.06); %m/s  %based on relationship in Liu et al. 2022 PNAS ('log' is natural logarithm)

%choice of k empirical model
kModel=questdlg('Which k model?','Gas transfer velocity','Raymond','Ulseth','Raymond');


%% 2. Preliminary calculations (Sc, k, atmospheric gas)

% start loop with 3 gases and assign concentration and flux variables
gases={'CO2','CH4','N2O'};

for g=1:3
    gas=gases{g};
    disp(' ');
    disp(['starting analysis for ', gas,'...'])
    disp(' ');
    switch gas
        case 'CO2'
            conc_field='CO2_converted';
            flux_field='F_CO2_converted';
            threshold=10;
        case 'CH4'
            conc_field='CH4_converted';
            flux_field='F_CH4_diff_converted';
            threshold=10;
        case 'N2O'
            conc_field='N2O_converted';
            flux_field='F_N2O_converted';
            threshold=5;
    end

    % calculate Schmidt numbers
    %first we gap fill missing temp data with mean annual air temp
    T.Temp_gapfill=T.Temperature;
    b=isnan(T.Temp_gapfill); T.Temp_gapfill(b)=T.airTemperature_BasAt(b);
    %then we get all 3 Sc numbers
    T.Sc_CO2=1923.6-125.06.*T.Temp_gapfill+4.3773.*T.Temp_gapfill.^2-0.085681.*T.Temp_gapfill.^3+0.00070284.*T.Temp_gapfill.^4; %Jahne et al (1987) / Wanninkhof (2014)
    T.Sc_CH4=1909.4-120.78.*T.Temp_gapfill+4.1555.*T.Temp_gapfill.^2-0.080578.*T.Temp_gapfill.^3+0.00065777.*T.Temp_gapfill.^4; %Jahne et al (1987) / Wanninkhof (2014)
    T.Sc_N2O=2141.2-152.56.*T.Temp_gapfill+5.8963.*T.Temp_gapfill.^2-0.12411.*T.Temp_gapfill.^3+0.0010655.*T.Temp_gapfill.^4; %Jahne et al (1987) / Wanninkhof (2014)

    % import coefficients of power-law relationship between ebullitive and
    % diffusive CH4 fluxes to add ebullitive component to all diffusive estimates
    if strcmp(gas,'CH4')
        p1=1.0636; p2=10^0.1850;
    end

    % calculate energy dissipation rate eD and k600 based on Ulseth et al. 2019 or Raymond et al. 2012
    T.eD=gravit.*T.slope_GRADES.*T.velocity; %eD in m2/s3
    switch kModel
        case 'Raymond'
            for a=1:size(T.eD,1)
                T.k600_modelled(a)=T.slope_GRADES(a).*T.velocity(a)*2841+2.02;
            end
        case 'Ulseth'
            for a=1:size(T.eD,1)
                T.k600_modelled(a)=T.slope_GRADES(a).*T.velocity(a)*2841+2.02;
                if T.eD(a)<0.02
                    T.k600_modelled(a)=exp(0.35*log(T.eD(a))+3.1); %k600 in m/d
                else
                    T.k600_modelled(a)=exp(1.18*log(T.eD(a))+6.43); %k600 in m/d
                end
            end
    end
    T.k600_modelled(T.k600_modelled==0)=0.5; %assign a minimum value for k600

    T=movevars(T,"k600_modelled",'After',"k600_converted");
    T=movevars(T,"slope_H90",'After',"Slope");
    T=movevars(T,"slope_GRADES",'After',"slope_H90");

    % calculate kCO2, kCH4 and kN2O
    T.(['k' gas])=T.k600_modelled.*(T.(['Sc_' gas])/600).^-0.5;

    % calculate atmospheric eq
    T.Temp_Kelvin=273.15+T.Temp_gapfill;
    T.Kh_CO2=10.^-(-(108.3865+0.01985076.*(T.Temp_Kelvin)-6919.53./(T.Temp_Kelvin)-40.4515.*log10(T.Temp_Kelvin)+669365./(T.Temp_Kelvin).^2)); %Plummer & Busenberg 1982 (10.1016/0016-7037(82)90056-4)
    T.Pressure(~isnan(T.elevation_H90))=(1-(.0000225577.*T.elevation_H90(~isnan(T.elevation_H90)))).^5.25588;
    T.Kh_CH4=exp(-68.8862+101.4956.*100./(T.Temp_Kelvin)+28.7314.*log((T.Temp_Kelvin)./100))./((T.Temp_Kelvin).*gas_constant./T.Pressure./100); %Wiesenburg & Guinasso 1979 (10.1021/je60083a006) (taken from Wanninkhof 2014)
    T.Kh_N2O=exp(-62.7062+97.3066.*100./(T.Temp_Kelvin)+24.1406.*log((T.Temp_Kelvin)./100))./((T.Temp_Kelvin).*gas_constant./T.Pressure./100); %Weiss & Price 1980 (10.1016/0304-4203(80)90024-9) (taken from Wanninkhof 2014)
    T.CO2_atm=CO2_sat.*T.Kh_CO2;
    T.CH4_atm=CH4_sat.*T.Kh_CH4;
    T.N2O_atm=N2O_sat.*T.Kh_N2O;



    %% 3. Estimate fluxes based on concentration and k

    % model fluxes when concentration data are available
    F=NaN(size(T,1),1); %preallocate
    for j=1:size(T,1)
        if ~isnan(T.(conc_field)(j))
            F(j)=(T.(conc_field)(j)-T.([gas '_atm'])(j)).*T.(['k' gas])(j); %F in mmol/m2/d
        else
            F(j)=NaN;
        end
    end
    T.f_modelled=F; clear F

    % check goodness of fit
    B=T(~isnan(T.f_modelled) & ~isnan(T.(flux_field)),:);
    B=B(B.f_modelled>0 & B.(flux_field)>0,:);
    [~,gof]=fit(B.f_modelled,B.(flux_field),'power1');
    disp('goodness of fit for measured/modelled flux data')
    disp(gof)
    disp(' ');
    clear B

    % create new column that merges modelled and measured f data
    for j=1:size(T,1)
        if ~isnan(T.(flux_field)(j))
            T.f_measmod(j)=T.(flux_field)(j);
        else
            T.f_measmod(j)=T.f_modelled(j);
        end
    end



    %% 4. Bootstrapping to obtain median fluxes + IQR per stream order and climate zone

    % create reduced datasets for the stream order / climate zone of interest
    f_measmod=NaN(boot_nb,9,5); sizes=NaN(9,5); %preallocate
    f_measmod_med_iqr=NaN(2,9,5); %preallocate
    for so=0:8 %this is the loop for stream orders
        for cl=1:5 %this is the loop for climate zones
            if so>=0 && so<8
                Sub=T(double(T.KG)==cl & T.shifted_streamOrder_H90==so,:); %create reduced dataset
                Sub_tot=T(T.shifted_streamOrder_H90==so,:); %entire dataset in case we don't have enough data in this climate zone
            elseif so==8
                Sub=T(double(T.KG)==cl & T.shifted_streamOrder_H90>=so,:); %the last category is all SO>9
                Sub_tot=T(T.shifted_streamOrder_H90>=so,:); %entire dataset in case we don't have enough data in this climate zone
            end

            %extract all fluxes, and save flux array for future use
            f_d_measmod=Sub.f_measmod(~isnan(Sub.f_measmod));
            f_d_measmod_tot=Sub_tot.f_measmod(~isnan(Sub_tot.f_measmod));
            if size(f_d_measmod,1)>=10
                if strcmp(gas,'CH4')
                    d=f_d_measmod+(p2.*f_d_measmod.^p1);
                    medFlux{so+1,cl}=d;
                    clear d
                else
                    medFlux{so+1,cl}=f_d_measmod;
                end
            else
                medFlux{so+1,cl}=NaN;
            end

            %bootstrap
            if size(f_d_measmod,1)>=threshold %this is the threshold above which we use the entire dataset - lower for N2O as there are less data for this gas. 
                data=f_d_measmod;
            else
                data=f_d_measmod_tot; %if not enough values we calculate a median across all climates
            end
            sizes(so+1,cl)=size(data,1);
            f_measmod(:,so+1,cl)=bootstrp(boot_nb,@median,data);
            f_measmod_med_iqr(1:2,so+1,cl)=bootci(boot_nb,{@median,data},'Alpha',0.5,'Type','bca'); %50% CI is equivalent to Q1 and Q3
            f_measmod_med_iqr(3,so+1,cl)=f_measmod_med_iqr(2,so+1,cl); %upper bound moved to 3rd row
            f_measmod_med_iqr(2,so+1,cl)=median(f_measmod(:,so+1,cl)); %median
            clear Sub f_d_measmod f_d_measmod_tot data
        end
    end

    % add ebullitive CH4 flux estimates based on power-law
    if strcmp(gas,'CH4')
        f_measmod_med_iqr=f_measmod_med_iqr+(p2.*f_measmod_med_iqr.^p1); %adding to bootstrapped values
        for j=1:size(T,1)
            if T.f_measmod(j)>0
                T.f_measmod(j)=T.f_measmod(j)+(p2.*T.f_measmod(j).^p1); %adding to R table
            end
        end
    end

    % output median fluxes
    disp(' ');
    disp(['median ',gas,' fluxes in mmol/m2/d across stream orders 1 (left) to 9+ (right)'])
    result='';
    for cl=1:5
        disp(climate_class{cl})
        for so=1:9
            val=medFlux{so,cl};
            md=median(val);
            switch gas
                case 'CO2'
                    result=sprintf('%.0f\t',md);
                case 'CH4'
                    result=sprintf('%.1f\t',md);
                case 'N2O'
                    result=sprintf('%.2f\t',md);
            end
            fprintf('%s',result);
            clear val
        end
        fprintf('\n');
        pause(1)
    end
    disp(' ')
    clear medFlux



    %% 5. Upscaling using river surface areas

    % estimate upscaled fluxes (F)
    F=NaN(boot_nb,9,5); F_med_iqr=NaN(3,9,5);
    for so=0:8
        for cl=1:5
            B=SurfAr(SurfAr.order==so & SurfAr.climate_class==cl,:); %find matching row in surface area dataset
            if isempty(B)==0 && ~isnan(f_measmod_med_iqr(2,so+1,cl))
                area_m2=B.surface_area_km2*1E6*(1-B.fraction_dry); %m2
                area_upper_m2=B.surface_area_km2*(1+percent_error_area/100)*1E6*(1-B.fraction_dry); % upper bound
                area_lower_m2=B.surface_area_km2*(1-percent_error_area/100)*1E6*(1-B.fraction_dry); %lower bound
                F(:,so+1,cl)=f_measmod(:,so+1,cl)*area_m2; %mmol/d
                F_med_iqr(2,so+1,cl)=f_measmod_med_iqr(2,so+1,cl)*area_m2; %median (mmol/d)
                F_med_iqr(1,so+1,cl)=f_measmod_med_iqr(1,so+1,cl)*area_lower_m2; %lower uncertainty bound (mmol/d)
                F_med_iqr(3,so+1,cl)=f_measmod_med_iqr(3,so+1,cl)*area_upper_m2; %upper uncertainty bound (mmol/d)
                clear B area_m2 area_upper_m2 area_lower_m2
            else
                F(:,so+1,cl)=NaN;
                F_med_iqr(:,so+1,cl)=NaN;
            end
        end
    end

    %convert fluxes to Tg C/yr and Tg N/yr 
    switch gas
        case {'CO2','CH4'}
            F=F.*molar_mass_C/1E15; % 1E15 mg/d = 1 Tg/d
            F=F*365; %Tg C /yr
            F_med_iqr=F_med_iqr.*molar_mass_C/1E15;
            F_med_iqr=F_med_iqr*365; %Tg C /yr
        case 'N2O'
            F=F.*molar_mass_N/1E15; % 1E15 mg/d = 1 Tg/d
            F=F*365; %Tg N /yr
            F_med_iqr=F_med_iqr.*molar_mass_N/1E15;
            F_med_iqr=F_med_iqr*365; %Tg N /yr
    end

    % calculate warming potentials
    switch gas
        case 'CO2'
            F_med_iqr_eq=F_med_iqr*molar_mass_CO2/molar_mass_C; % Tg CO2 yr-1
        case 'CH4'
            F_med_iqr_gas=F_med_iqr*molar_mass_CH4/molar_mass_C; % Tg CH4 yr-1
            F_med_iqr_eq=F_med_iqr_gas*GWP_CH4; % Tg CO2-eq yr-1
        case 'N2O'
            F_med_iqr_gas=F_med_iqr*molar_mass_N2O/molar_mass_N; % Tg N2O yr-1
            F_med_iqr_eq=F_med_iqr_gas*GWP_N2O; % Tg CO2-eq yr-1
    end

            % output total median flux for a given gas
        for i=1:3 %lower, med, upper
            for cl=1:5
                sums_matrix(i,cl)=sum(F_med_iqr(i,:,cl),'omitnan'); %Tg C/yr or Tg N/yr
            end
        end
        switch gas
            case 'CO2'
                sums_matrix=sums_matrix*molar_mass_CO2/molar_mass_C; %Tg CO2/yr
            case 'CH4'
                sums_matrix=sums_matrix*molar_mass_CH4/molar_mass_C; %Tg CH4/yr
            case 'N2O'
                sums_matrix=sums_matrix*molar_mass_N2O/molar_mass_N; %Tg N2O/yr
        end
        disp(' ')
        disp(['total median ',gas,' fluxes in Tg gas/yr across climate classes'])
        disp('(last line is total)')
        for col=1:size(sums_matrix,2)
            dn=sums_matrix(1,col);
            md=sums_matrix(2,col);
            up=sums_matrix(3,col);
            switch gas
                case 'CO2'
                    result=sprintf('%.0f (%.0f-%.0f)',md,dn,up);
                case 'CH4'
                    result=sprintf('%.2f (%.2f-%.2f)',md,dn,up);
                case 'N2O'
                    result=sprintf('%.3f (%.3f-%.3f)',md,dn,up);
            end
            disp(result)
        end
        total_dn=sum(sums_matrix(1,:));
        total_md=sum(sums_matrix(2,:));
        total_up=sum(sums_matrix(3,:));
        switch gas
            case 'CO2'
                total_result=sprintf('%.0f (%.0f-%.0f)',total_md,total_dn,total_up);
            case 'CH4'
                total_result=sprintf('%.2f (%.2f-%.2f)',total_md,total_dn,total_up);
            case 'N2O'
                total_result=sprintf('%.3f (%.3f-%.3f)',total_md,total_dn,total_up);
        end
        disp(total_result)
        disp (' ')
        pause(2)


    %save data before starting another gas
    arealFlux_=['arealFlux_' gas];
    scaledFlux_=['scaledFlux_' gas];
    scaledFlux_eq_=['scaledFlux_eq_' gas];

    eval([arealFlux_ ' = f_measmod_med_iqr;']);
    eval([scaledFlux_ ' = F_med_iqr;']);
    eval([scaledFlux_eq_ ' = F_med_iqr_eq;']);

end


disp(' ')
disp('flux calculation for streams/rivers done')
disp(' ')

end


