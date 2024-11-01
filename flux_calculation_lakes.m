%% Function to upscale GHG fluxes per lake size class and climate zone using a bootstrapping approach
%
%
function[scaledFlux_eq_CO2_lake,scaledFlux_eq_CO2_reservoir,...
    scaledFlux_eq_CH4_lake,scaledFlux_eq_CH4_reservoir,...
    scaledFlux_eq_N2O_lake,scaledFlux_eq_N2O_reservoir]=flux_calculation_lakes(T)
%
% This function contains the following sections:
% 1. Upload data and constants
% 2. Flux calculation via bootstrapping
% 3. Upscaling using lake/reservoir surface areas
%
% Inputs
% 'T_per_system' is a reduced table that has been processed using 'unit_conversion', 'climate_classification', 'medians_per_site' and 'weighted means per lake'
% In addition, the following files are needed in the input dataset folder:
% 'grouped_lakes_by_climate_and_area_natural.csv'
% 'grouped_lakes_by_climate_and_area_reservoirs.csv'
%
% Outputs
% scaledFlux_eq_gas_system is a 3x6x5 array with Rows: Q1, med, Q3; Columns: size classes; Depths: climate zones.
% Values in scaledFlux_eq_gas_system are upscaled fluxes in Tg CO2-eq/yr.



%% 1. Upload data and constants

disp(' '); disp('starting lake/reservoir flux calculation...'); disp(' ')

%constants
boot_nb=1000; %nb of bootstrap samples
molar_mass_C=12.011; %g mol-1
molar_mass_N=14.01; %g mol-1
molar_mass_N2O=44.02; %g mol-1
molar_mass_CO2=44.01; %g mol-1
molar_mass_CH4=16.05; %g mol-1
GWP_N2O=298;
GWP_CH4=34;
percent_error_area=10; %percent error in surface area (10% for standing waters)

%import lake and reservoir surface areas for each size and climate class
A1=readtable('grouped_lakes_by_climate_and_area_natural.csv');
A1.new_climate_class=NaN(height(A1),1);
A1.new_climate_class(ismember(A1.climate_zone,[1,2]))=1;
A1.new_climate_class(ismember(A1.climate_zone,3))=2;
A1.new_climate_class(ismember(A1.climate_zone,[4,5,6,7]))=3;
A1.new_climate_class(ismember(A1.climate_zone,[11,14]))=4;
A1.new_climate_class(ismember(A1.climate_zone,[12,15,22,23,29,30]))=5;
summary_table=groupsummary(A1,{'new_climate_class','area_class'},'sum','area_km2');
SurfArLakes=removevars(summary_table,'GroupCount');
SurfArLakes=SurfArLakes(~isnan(SurfArLakes.new_climate_class),:);
SurfArLakes=renamevars(SurfArLakes,'new_climate_class','climate_class');
clear summary_table
A2=readtable('grouped_lakes_by_climate_and_area_reservoirs.csv');
A2.new_climate_class=NaN(height(A2),1);
A2.new_climate_class(ismember(A2.climate_zone,[1,2]))=1;
A2.new_climate_class(ismember(A2.climate_zone,3))=2;
A2.new_climate_class(ismember(A2.climate_zone,[4,5,6,7]))=3;
A2.new_climate_class(ismember(A2.climate_zone,[11,14]))=4;
A2.new_climate_class(ismember(A2.climate_zone,[12,15,22,23,29,30]))=5;
summary_table=groupsummary(A2,{'new_climate_class','area_class'},'sum','area_km2');
SurfArRes=removevars(summary_table,'GroupCount');
SurfArRes=SurfArRes(~isnan(SurfArRes.new_climate_class),:);
SurfArRes=renamevars(SurfArRes,'new_climate_class','climate_class');
clear summary_table
SurfAr={SurfArLakes;SurfArRes};

% import coefficients of power-law relationship between ebullitive and diffusive CH4 fluxes
p1=1.0636; p2=10^0.1850;



%% 2. Flux calculation via bootstrapping

% separate between reservoirs and lakes/ponds
mLake=T(strcmp(T.type,'lake') | strcmp(T.type,'pond'),:);
mReservoir=T(strcmp(T.type,'reservoir'),:);
Dbs={mLake; mReservoir};
type={'lake','reservoir'};

% start loop with 3 gases and assign flux variables
gases={'CO2','CH4','N2O'};

for g=1:3
    gas=gases{g};
    disp(' ');
    disp(['starting analysis for ', gas,'...'])
    disp(' ');
    switch gas
        case 'CO2'
            flux_field='weighted_FCO2';
        case 'CH4'
            flux_field='weighted_FCH4_diff';
        case 'N2O'
            flux_field='weighted_FN2O';
    end

    % create reduced datasets for the system type / size class / climate zone of interest
    f=NaN(boot_nb,5,5);  f_ci=NaN(2,5,5); nb=NaN(1,5,5); %preallocate
    for resnat=1:2 %this is the loop for lakes/reservoirs
        M=Dbs{resnat}; %first we run lakes, then reservoirs
        for lks=1:6 %this is the loop for lake size classes
            for cl=1:5 %this is the loop for climate zones
                if lks==1
                    S=M(M.lake_area>0.01 & M.lake_area<0.1,:);
                elseif lks==2
                    S=M(M.lake_area>0.1 & M.lake_area<1,:);
                elseif lks==3
                    S=M(M.lake_area>1 & M.lake_area<10,:);
                elseif lks==4
                    S=M(M.lake_area>10 & M.lake_area<100,:);
                elseif lks==5
                    S=M(M.lake_area>100 & M.lake_area<1000,:);
                elseif lks==6 && resnat==1
                    S=M(M.lake_area>1000,:);
                elseif lks==6 && resnat==2
                    S=M(M.lake_area>100 & M.lake_area<1000,:); %special condition for reservoirs >1000: we use data from reservoirs 100-1000 (see Methods).
                end

                %get flux data and remove NaNs
                f_d_tot=S.(flux_field)(~isnan(S.(flux_field)));

                %bootstrapping
                if size(f_d_tot,1)>1
                    f(:,lks,cl)=bootstrp(boot_nb,@median,f_d_tot);
                    f_ci(:,lks,cl)=bootci(boot_nb,{@median,f_d_tot},'Alpha',0.5,'Type','bca'); %50% CI is equivalent to Q1 and Q3
                    nb(:,lks,cl)=size(f_d_tot,1);
                else
                    f(:,lks,cl)=NaN;
                    f_ci(:,lks,cl)=NaN;
                    nb(:,lks,cl)=size(f_d_tot,1);
                end

                clear Sub f_d_tot S
            end
        end

        % add ebullitive CH4 flux estimates based on power-law
        if strcmp(gas,'CH4')
            f_ci=f_ci+(p2.*f_ci.^p1);
            f=f+(p2.*f.^p1);
        end


        
        %% 3. Upscaling using lake/reservoir surface areas
        % multiplication of areal fluxes by surface area
        F=NaN(boot_nb,6,5); F_ci=NaN(3,6,5); f_meds=NaN(6,5); %preallocate
        Ql=SurfAr{resnat}; %this is the surface area data
        for lks=1:6
            for cl=1:5
                B=Ql(Ql.area_class==lks+1 & Ql.climate_class==cl,:); %find matching row
                if lks==6
                    B2=Ql(Ql.area_class==lks+2 & Ql.climate_class==cl,:); %for last size class we add >10000 km2
                    B.sum_area_km2=B.sum_area_km2+B2.sum_area_km2;
                end
                area_m2=B.sum_area_km2*1E6; %m2
                area_upper_m2=B.sum_area_km2*(1+0.01*percent_error_area)*1E6; % upper bound
                area_lower_m2=B.sum_area_km2*(1-0.01*percent_error_area)*1E6; %lower bound

                F(:,lks,cl)=f(:,lks,cl)*area_m2; %mmol/d
                F_ci(2,lks,cl)=median(F(:,lks,cl)); %mmol/d
                F_ci(1,lks,cl)=f_ci(1,lks,cl)*area_lower_m2; %lower uncertainty bound (mmol/d)
                F_ci(3,lks,cl)=f_ci(2,lks,cl)*area_upper_m2; %upper uncertainty bound (mmol/d)

                f_meds(lks,cl)=F_ci(2,lks,cl)/B.sum_area_km2/1E6; %mmol/m2/d
                clear B B2
            end
        end
        clear Ql

        % output median fluxes for paper
        disp(['median ',gas,' flux in mmol/m2/d for ',type{resnat},' systems across size classes'])
        disp(f_meds(:,1))
        pause(2)

        %convert fluxes to Tg C/yr and Tg N/yr
        switch gas
            case {'CO2','CH4'}
                F=F.*molar_mass_C/1E15; % 1E15 mg/d = 1 Tg/d
                F=F*365; %Tg C /yr
                F_ci=F_ci.*molar_mass_C/1E15;
                F_ci=F_ci*365; %Tg C /yr
            case 'N2O'
                F=F.*molar_mass_N/1E15; % 1E15 mg/d = 1 Tg/d
                F=F*365; %Tg N /yr
                F_ci=F_ci.*molar_mass_N/1E15;
                F_ci=F_ci*365; %Tg N /yr
        end

        % calculate warming potentials
        switch gas
            case 'CO2'
                F_ci_eq=F_ci*molar_mass_CO2/molar_mass_C; % Tg CO2 yr-1
            case 'CH4'
                F_ci_eq_gas=F_ci*molar_mass_CH4/molar_mass_C; % Tg CH4 yr-1
                F_ci_eq=F_ci_eq_gas*GWP_CH4; % Tg CO2-eq yr-1
            case 'N2O'
                F_ci_eq_gas=F_ci*molar_mass_N2O/molar_mass_N; % Tg N2O yr-1
                F_ci_eq=F_ci_eq_gas*GWP_N2O; % Tg CO2-eq yr-1
        end

        % output total median flux for a given gas
        for i=1:3 %lower, med, upper
            for cl=1:5
                sums_matrix(i,cl)=nansum(F_ci(i,:,cl)); %Tg C/yr or Tg N/yr
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
        disp(['total median ',gas,' fluxes in Tg gas/yr for ',type{resnat},' systems across climate classes'])
        disp('(last line is total)')
        for col=1:size(sums_matrix,2)
            dn=sums_matrix(1,col);
            md=sums_matrix(2,col);
            up=sums_matrix(3,col);
            switch gas
                case 'CO2'
                    result=sprintf('%.1f (%.1f-%.1f)',md,dn,up);
                case 'CH4'
                    result=sprintf('%.2f (%.2f-%.2f)',md,dn,up);
                case 'N2O'
                    result=sprintf('%.3f (%.3f-%.3f)',md,dn,up);
            end
            disp(result)
        end
        total_dn(resnat)=sum(sums_matrix(1,:));
        total_md(resnat)=sum(sums_matrix(2,:));
        total_up(resnat)=sum(sums_matrix(3,:));
        switch gas
            case 'CO2'
                total_result=sprintf('%.1f (%.1f-%.1f)',total_md(resnat),total_dn(resnat),total_up(resnat));
            case 'CH4'
                total_result=sprintf('%.2f (%.2f-%.2f)',total_md(resnat),total_dn(resnat),total_up(resnat));
            case 'N2O'
                total_result=sprintf('%.3f (%.3f-%.3f)',total_md(resnat),total_dn(resnat),total_up(resnat));
        end
        disp(total_result)
        disp (' ')
        pause(2)

        if resnat==2
            total_dn=sum(total_dn);
            total_md=sum(total_md);
            total_up=sum(total_up);

            disp(['Results Lakes + reservoirs for ',gas])
            switch gas
                case 'CO2'
                    total_result=sprintf('%.1f (%.1f-%.1f)',total_md,total_dn,total_up);
                case 'CH4'
                    total_result=sprintf('%.2f (%.2f-%.2f)',total_md,total_dn,total_up);
                case 'N2O'
                    total_result=sprintf('%.3f (%.3f-%.3f)',total_md,total_dn,total_up);
            end
            disp(total_result)
            pause(2)
            disp(' ')
        end

    %save data before starting another system type and gas
    scaledFlux_eq_=['scaledFlux_eq_' gas '_' type{resnat}];
    eval([scaledFlux_eq_ ' = F_ci_eq;']);

    end

end

disp(' ')
disp('flux calculation for lakes/reservoirs done')
disp(' ')

end





