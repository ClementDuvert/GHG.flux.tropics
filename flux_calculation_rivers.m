%% Function to upscale GHG fluxes by stream order and climate zone using a bootstrapping approach
% The function uses either the Raymond (eq5) or Ulseth gas transfer velocity model (depending on selection in the prompts). 
% The function computes either diffusive, ebullitive, or total CH4 fluxes (depending on selection in the prompts).
%
function[T,errorBar_CO2,errorBar_CH4,errorBar_N2O,arealFlux_CO2,scaledFlux_CO2,scaledFlux_eq_CO2,arealFlux_CH4,...
    scaledFlux_CH4,scaledFlux_eq_CH4,arealFlux_N2O,scaledFlux_N2O,...
    scaledFlux_eq_N2O]=flux_calculation_rivers(T)
%
% This function contains the following sections:
% 1. Initialise constants and import dataset
% 2. Preliminary calculations (Sc, k, atmospheric gas)
% 3. Estimate fluxes based on concentration and k
% 4. Bootstrapping to obtain median (Q1-Q3) fluxes per stream order and climate zone
% 5. Upscaling using river surface areas
%
% Inputs
% T is a reduced table that has been processed using 'unit_conversion', 'climate_classification' and 'averaging_per_site'
% In addition, the following file is needed in the input dataset folder:
% 'river_surface_areas.csv'
%
% Outputs
% all outputs are organised as follows. 3 rows: Q1, med, Q3. 9 columns: stream order 1 to 9+. 5 depths: climate zones.
% arealFlux_gas is a 3x9x5 array with areal fluxes in mmol/m2/d.
% scaledFlux_gas is is a 3x9x5 array with upscaled fluxes in Tg C/yr or Tg N/yr.
% scaledFlux_eq_gas is a 3x9x5 array with upscaled fluxes in Tg CO2-eq/yr.


%% 1. Initialise constants and import surface area dataset

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
threshold=10;

%import river surface areas for each SO and climate class
SurfAr=readtable('river_surface_areas.csv');

%choice of k empirical model
kModel=questdlg('Which k model?','Gas transfer velocity','Raymond','Ulseth','Raymond');

%choice of mean/median for bootstrapping
meanmed=questdlg('Use mean or median for bootstrapping?','bootstrapping','mean','median','median');

%choice of diffusive or ebullitive or total flux for CH4
methane=questdlg('which methane flux?','methane','diffusive','ebullitive','total','diffusive');

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
        case 'CH4'
            conc_field='CH4_converted';
            flux_field='F_CH4_diff_converted';
        case 'N2O'
            conc_field='N2O_converted';
            flux_field='F_N2O_converted';
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
        %p1=1.0636; p2=10^0.1850;
        p1=0.9123; p2=10^0.2658;   %revised May 2025
    end

    % calculate energy dissipation rate eD and k600 based on Ulseth et al. 2019 or Raymond et al. 2012
    T.eD=gravit.*T.slope_consolidated.*T.velocity; %eD in m2/s3
    switch kModel
        case 'Raymond'
            for a=1:size(T.eD,1)
                T.k600_modelled(a)=T.slope_consolidated(a).*T.velocity(a)*2841+2.02; %k600 in m/d
            end
        case 'Ulseth'
            for a=1:size(T.eD,1)
                if T.eD(a)<0.02
                    T.k600_modelled(a)=exp(0.35*log(T.eD(a))+3.1); %k600 in m/d
                else
                    T.k600_modelled(a)=exp(1.18*log(T.eD(a))+6.43); %k600 in m/d
                end
            end
    end

    T.k600_modelled(T.k600_modelled==0)=0.5; %assign a minimum value for k600
    T=movevars(T,"k600_modelled",'After',"k600_converted");
    T.k600_consolidated=T.k600_converted; %we always prefer reported k600
    T.k600_consolidated(isnan(T.k600_converted))=T.k600_modelled(isnan(T.k600_converted)); %if no reported k600, use modelled one
    T=movevars(T,"k600_consolidated",'After',"k600_modelled");

    % calculate kCO2, kCH4 and kN2O
    T.(['k' gas])=T.k600_consolidated.*(T.(['Sc_' gas])/600).^-0.5;

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

    % check goodness of fit for measured vs modelled flux
    B2=T(~isnan(T.f_modelled) & ~isnan(T.(flux_field)) & T.F_method_cat=='measured',:); %we only keep 'truly' measured fluxes
    B2=B2(B2.f_modelled>0 & B2.(flux_field)>0,:);
    GOF_mod{g}=B2.f_modelled;
    GOF_meas{g}=B2.(flux_field);

    if g==3
        all_GOF_mod=vertcat(GOF_mod{:});
        all_GOF_meas=vertcat(GOF_meas{:});
        log_mod=log10(all_GOF_mod);
        log_meas=log10(all_GOF_meas);
        n1=numel(all_GOF_meas);
        curve=fitlm(log_meas,log_mod,'poly1');
        p_val=curve.Coefficients.pValue(2);
        coeffs=curve.Coefficients.Estimate;
        intercept=coeffs(1);
        slope=coeffs(2);
        predicted=predict(curve,log_meas);
        residuals=log_mod-predicted;
        rmse_orig=sqrt(mean((10.^residuals).^2));
      
disp('goodness of fit for measured/modelled flux data')
        disp(curve.Coefficients)
        disp(' ');
        figure('Position',[440 325 566 472]);

        subplot(2,2,1)
        xlim([10^-2.5;10^4]); ylim([10^-2.5;10^4])
        xlims=xlim; ylims=ylim;
        plot([min(xlims(1),ylims(1)),max(xlims(2),ylims(2))],[min(xlims(1),ylims(1)),max(xlims(2),ylims(2))],'k--','LineWidth',1.2);
        hold on
        scatter(all_GOF_meas,all_GOF_mod,60,[0.2 0.2 0.2],'filled','MarkerFaceAlpha',0.3);
        set(gca,'XScale','log','YScale','log')
        x_fit=[min(all_GOF_meas) max(all_GOF_meas)];
        y_fit=10.^(slope*log10(x_fit)+intercept);
        plot(x_fit,y_fit,'-r','Color',[.86 .2 .2],'LineWidth',1);
        xlabel('measured flux (mmol m^{-2} d^{-1})');
        ylabel('modelled flux (mmol m^{-2} d^{-1})');
        set(gca,'XTick',[1e-2 1e0 1e2 1e4],'XTickLabel',{'10^{-2}','10^{0}','10^{2}','10^{4}'})
        set(gca,'YTick',[1e-2 1e0 1e2 1e4],'YTickLabel',{'10^{-2}','10^{0}','10^{2}','10^{4}'})
        grid on; box on
        set(gca,'FontSize',12,'FontName','Arial');
        a=10^intercept;  %coefficient in original space
        b=slope;    %exponent
        eqn_txt=sprintf('f_{mod}=%.2f × f_{meas}^{%.2f}',a,b);
        stat_txt=sprintf('R^2=%.2f; RMSE=%.1f',curve.Rsquared.Adjusted,rmse_orig);
        n_txt=sprintf('n=%d; p=%.1e',n1,p_val);
        text(0.05,0.9,stat_txt,'Units','normalized','FontSize',9,'Color',[.86 .2 .2]);
        text(0.05,0.82,eqn_txt,'Units','normalized','FontSize',9,'Color',[.86 .2 .2]);
        text(0.05,0.74,n_txt,'Units','normalized','FontSize', 9,'Color',[.86 .2 .2]);
        xlim([10^-2.5;10^4]); ylim([10^-2.5;10^4])

        subplot(2,2,3)
        scatter(all_GOF_meas,residuals,60,[0.2 0.2 0.2],'filled','MarkerFaceAlpha',0.3);
        hold on
        set(gca,'XScale','log')
        set(gca,'XTick',[1e-2 1e0 1e2 1e4],'XTickLabel',{'10^{-2}','10^{0}','10^{2}','10^{4}'})
        xlabel('measured flux (mmol m^{-2} d^{-1})');
        ylabel('log(residuals)');
        plot([1e-5 1e4],[0 0],'k--','LineWidth',1.2);
        grid on; box on
        set(gca,'FontSize',12,'FontName','Arial');
        xlim([10^-2.5;10^4])
        clear B2
    end

    % create new column that merges modelled and measured f data
    measmod_field=['f_measmod_' gas];
    T.(measmod_field)=NaN(size(T,1),1);
    for j=1:size(T,1)
        if ~isnan(T.(flux_field)(j))
            T.(measmod_field)(j)=T.(flux_field)(j);
        else
            T.(measmod_field)(j)=T.f_modelled(j);
        end
    end


    %% 4. Bootstrapping to obtain median fluxes + IQR per stream order and climate zone

    % create reduced datasets for the stream order / climate zone of interest
    f_measmod=NaN(boot_nb,9,5); sizes=NaN(9,5); %preallocate
    f_measmod_med_iqr=NaN(2,9,5); %preallocate
    for so=1:9 %this is the loop for stream orders
        for cl=1:5 %this is the loop for climate zones
            if so>=1 && so<9
                Sub=T(double(T.KG)==cl & T.streamorder_consolidated==so,:); %create reduced dataset
                Sub_tot=T(T.streamorder_consolidated==so,:); %entire dataset in case we don't have enough data in this climate zone
            elseif so==9
                Sub=T(double(T.KG)==cl & T.streamorder_consolidated>=so,:); %the last category is all SO>9
                Sub_tot=T(T.streamorder_consolidated>=so,:); %entire dataset in case we don't have enough data in this climate zone
            end

            %extract all fluxes, and save flux array for future use
            f_d_measmod=Sub.(measmod_field)(~isnan(Sub.(measmod_field)));
            f_d_measmod_tot=Sub_tot.(measmod_field)(~isnan(Sub_tot.(measmod_field)));
            medFlux{so,cl}=f_d_measmod;

            %bootstrap
            if size(f_d_measmod,1)>=threshold %this is the threshold above which we use the entire dataset
                data=f_d_measmod;
            else
                data=f_d_measmod_tot; %if not enough values we calculate a median across all climates
                disp(['all data used for subgroup climate ',num2str(cl),' and SO ',num2str(so)])
            end
            sizes(so,cl)=size(data,1);
             switch meanmed
                    case 'mean'
            f_measmod(:,so,cl)=bootstrp(boot_nb,@mean,data);
            f_measmod_med_iqr(1:2,so,cl)=bootci(boot_nb,{@mean,data},'Alpha',0.5,'Type','bca'); %50% CI is equivalent to Q1 and Q3
            f_measmod_med_iqr(3,so,cl)=f_measmod_med_iqr(2,so,cl); %upper bound moved to 3rd row
            f_measmod_med_iqr(2,so,cl)=mean(f_measmod(:,so,cl)); %mean
                 case 'median'
                         f_measmod(:,so,cl)=bootstrp(boot_nb,@median,data);
            f_measmod_med_iqr(1:2,so,cl)=bootci(boot_nb,{@median,data},'Alpha',0.5,'Type','bca'); %50% CI is equivalent to Q1 and Q3
            f_measmod_med_iqr(3,so,cl)=f_measmod_med_iqr(2,so,cl); %upper bound moved to 3rd row
            f_measmod_med_iqr(2,so,cl)=median(f_measmod(:,so,cl)); %median                    
             end             
                     clear Sub f_d_measmod f_d_measmod_tot data
        end
    end

    % add ebullitive CH4 flux estimates based on power-law if required
    if strcmp(gas,'CH4')
        switch methane
            case 'diffusive'
                %no transformation needed
            case 'ebullitive'
                f_measmod_med_iqr=p2.*f_measmod_med_iqr.^p1; %adding to bootstrapped values
                for j=1:size(T,1)
                    if T.(measmod_field)(j)>0
                        T.(measmod_field)(j)=p2.*T.(measmod_field)(j).^p1; %adding to table
                    end
                end
            case 'total'
                f_measmod_med_iqr=f_measmod_med_iqr+(p2.*f_measmod_med_iqr.^p1); %adding to bootstrapped values
                for j=1:size(T,1)
                    if T.(measmod_field)(j)>0
                        T.(measmod_field)(j)=T.(measmod_field)(j)+(p2.*T.(measmod_field)(j).^p1); %adding to table
                    end
                end
        end
    end

% Output formatted table: "median (IQR; n=X)" across stream orders
disp(' ')
disp(['Median ',gas,' fluxes (',meanmed, ') in mmol/m2/d across stream orders 1–9+'])
disp('Format: median (IQR; n=X)')
disp(' ')

for cl=1:5
    fprintf('%-12s\t',climate_class{cl});
    for so=1:9
        val=medFlux{so,cl};
        val=val(~isnan(val));
        n=numel(val);
        if n==0
            result='NA';
        else
            switch meanmed
                case 'mean'
                    central=mean(val);
                    q1=prctile(val,25);
                    q3=prctile(val,75);
                case 'median'
                    central=median(val);
                    q1=prctile(val,25);
                    q3=prctile(val,75);
            end
            switch gas
                case 'CO2'
                    result=sprintf('%.1f (%.1f–%.1f; n=%d)',central,q1,q3,n);
                case 'CH4'
                    result=sprintf('%.3f (%.3f–%.3f; n=%d)',central,q1,q3,n);
                case 'N2O'
                    result=sprintf('%.3f (%.3f–%.3f; n=%d)',central,q1,q3,n);
            end
        end

        fprintf('%-28s', result);  % print per stream order
        clear val
    end
    fprintf('\n');
end

disp(' ')
clear medFlux


    %% 5. Upscaling using river surface areas

    % estimate upscaled fluxes (F)
    F=NaN(boot_nb,9,5); F_med_iqr=NaN(3,9,5);
    for so=0:8 %in this table the stream orders 1 are assigned 0s, and so on.
        for cl=1:5
            B=SurfAr(SurfAr.order==so & SurfAr.climate_class==cl,:); %find matching row in surface area dataset
            if isempty(B)==0 && ~isnan(f_measmod_med_iqr(2,so+1,cl)) %the +1 is to match the right SO (1 for 0 in area file)
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

    % calculate error bars per climate and gas
    F_ci_eq_errorbar=NaN(5,2);
    for cl=1:5
        F_ci_eq_errorbar(cl,1)=sum(F_med_iqr_eq(1,:,cl),'omitnan');
        F_ci_eq_errorbar(cl,2)=sum(F_med_iqr_eq(3,:,cl),'omitnan');
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
    pause(1)

    %save data before starting another gas
    errorBar_=['errorBar_' gas];
    arealFlux_=['arealFlux_' gas];
    scaledFlux_=['scaledFlux_' gas];
    scaledFlux_eq_=['scaledFlux_eq_' gas];

    eval([errorBar_ ' = F_ci_eq_errorbar;']);
    eval([arealFlux_ ' = f_measmod_med_iqr;']);
    eval([scaledFlux_ ' = F_med_iqr;']);
    eval([scaledFlux_eq_ ' = F_med_iqr_eq;']);

    disp('subgroup sizes')
    disp(sizes)

end

disp(' ')
disp('flux calculation for streams/rivers done')
disp(' ')

end


