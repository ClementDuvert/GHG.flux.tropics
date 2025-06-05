%% Function to homogenise units for concentrations, fluxes, gas transfer velocities, and discharge
%
function[Tclean]=unit_conversion(T,choice)
%
% Inputs
%  - T must be a table with all variables listed as per txt file
%  - choice is 'streams/rivers' or 'lakes/reservoirs'


%% 1. Initialise constants and variables
disp(' '); disp('starting unit conversion...'); disp(' ')

% constants
gas_constant=8.3144626181; %m3 Pa mol-1 K-1
molar_mass_C=12.011; %g mol-1
molar_mass_N=14.01; %g mol-1
molar_mass_CO2=44.01; %g mol-1
molar_mass_CH4=16.04; %g mol-1
molar_mass_N2O=44.013; %g mol-1
molar_mass_O2=31.999; %g mol-1
CO2_sat=415; %umol/mol or ppmv  in 2020
CH4_sat=1.87; %umol/mol or ppmv  in 2020
N2O_sat=0.332; %umol/mol or ppmv  in 2020
hours_in_a_day=24;
minutes_in_a_day=1440;
seconds_in_a_day=86400;
days_in_a_year=365;
seconds_in_a_year=seconds_in_a_day*days_in_a_year;
seconds_in_a_minute=60;
days_in_a_month=30.437; %average
milligrams_in_a_tonne=10^9;
m2_in_a_ha=10000;
centimeters_in_a_meter=100;
litres_in_a_m3=1000;
m3_in_a_km3=1E9;
m3_in_a_ft3=0.0283168;
m3_in_a_ML=1000;
p0=1.01325; %standard atmospheric pressure at 0m in bar (=1 atm)

%pressure calculation
a=repelem(p0,size(T,1)); %assign pressure at elevation 0m for all data
T.Pressure=a';
if isequal(choice,'streams/rivers')
    T.Pressure(~isnan(T.elevation_H90))=(1-(.0000225577.*T.elevation_H90(~isnan(T.elevation_H90)))).^5.25588; %recalculate pressure based on H90 elevation (rivers)
elseif isequal(choice,'lakes/reservoirs')
    T.Pressure(~isnan(T.Elevation))=(1-(.0000225577.*T.Elevation(~isnan(T.Elevation)))).^5.25588; %recalculate pressure based on elevation (lakes)
end

% Henry constant calculation (in mol/L/atm)
T.Temp_gapfill=T.Temperature;
b=isnan(T.Temp_gapfill);
if isequal(choice,'streams/rivers')
    T.Temp_gapfill(b)=T.airTemperature_BasAt(b); %if water temperature is not reported, assign air temperature
elseif isequal(choice,'lakes/reservoirs')
        T.Temp_gapfill(b)=18; %for lakes we don't have air temperature; to be amended
end
totalFilled=sum(b);
disp(['gap-filled ',num2str(totalFilled),' temperature data out of ',num2str(size(T,1)),' total']);
T.Temp_Kelvin=273.15+T.Temp_gapfill;
%Plummer & Busenberg 1982 (10.1016/0016-7037(82)90056-4)
T.Kh_CO2=10.^-(-(108.3865+0.01985076.*(T.Temp_Kelvin)-6919.53./(T.Temp_Kelvin)-40.4515.*log10(T.Temp_Kelvin)+669365./(T.Temp_Kelvin).^2));
%Wiesenburg & Guinasso 1979 (10.1021/je60083a006) (taken from Wanninkhof 2014)
T.Kh_CH4=exp(-68.8862+101.4956.*100./(T.Temp_Kelvin)+28.7314.*log((T.Temp_Kelvin)./100))./((T.Temp_Kelvin).*gas_constant./T.Pressure./100);
%Weiss & Price 1980 (10.1016/0304-4203(80)90024-9) (taken from Wanninkhof 2014)
T.Kh_N2O=exp(-62.7062+97.3066.*100./(T.Temp_Kelvin)+24.1406.*log((T.Temp_Kelvin)./100))./((T.Temp_Kelvin).*gas_constant./T.Pressure./100);

% make the 'CO2 method' and 'Flux method' variables categorical
T=convertvars(T,'CO2_method_cat','categorical');
T=convertvars(T,'F_method_cat','categorical');


%% 2. Convert CO2 concentration units
%find unit type for each row in DB
units_uatm={'µatm','uatm','_atm','uatm at 25degC'}';
co2_1=cellfun(@(c)strcmp(c,units_uatm),T.CO2_units,'UniformOutput',false);
units_ppm={'ppm','ppmv'}';
co2_2=cellfun(@(c)strcmp(c,units_ppm),T.CO2_units,'UniformOutput',false);
units_umol={'_mol/L','_mol L-1','µM','_mol L_1','uM','_M','umol/L','µmol CO2/L','umolL-1'}';
co2_3=cellfun(@(c)strcmp(c,units_umol),T.CO2_units,'UniformOutput',false);
units_nmol={'nmol/L'}';
co2_4=cellfun(@(c)strcmp(c,units_nmol),T.CO2_units,'UniformOutput',false);
units_mmol={'mM','mmol/L'}';
co2_5=cellfun(@(c)strcmp(c,units_mmol),T.CO2_units,'UniformOutput',false);
units_mgL={'mg/L','mgC/L','mg C/L'}';
co2_6=cellfun(@(c)strcmp(c,units_mgL),T.CO2_units,'UniformOutput',false);
units_ugL={'ug/L','_g/L'}';
co2_7=cellfun(@(c)strcmp(c,units_ugL),T.CO2_units,'UniformOutput',false);
units_sat={'%sat','% sat','% saturation','%Sat'}';
co2_8=cellfun(@(c)strcmp(c,units_sat),T.CO2_units,'UniformOutput',false);
units_mgCO2L={'g CO2/m3','mg CO2/L'}';
co2_9=cellfun(@(c)strcmp(c,units_mgCO2L),T.CO2_units,'UniformOutput',false);
units_atm={'atm'}';
co2_10=cellfun(@(c)strcmp(c,units_atm),T.CO2_units,'UniformOutput',false);

%convert all CO2 units into umol/L
T.CO2_converted=NaN(size(T,1),1);
for x=1:size(T,1)
    if any(co2_1{x}) %uatm
        T.CO2_converted(x)=T.CO2(x)*T.Kh_CO2(x);
    elseif any(co2_2{x}) %ppm
        T.CO2_converted(x)=T.CO2(x)/T.Pressure(x)*T.Kh_CO2(x);
    elseif any(co2_3{x}) %umol/L
        T.CO2_converted(x)=T.CO2(x);
    elseif any(co2_4{x}) %nmol/L
        T.CO2_converted(x)=T.CO2(x)/1000;
    elseif any(co2_5{x}) %mmol/L
        T.CO2_converted(x)=T.CO2(x)*1000;
    elseif any(co2_6{x}) %mg/L
        T.CO2_converted(x)=T.CO2(x)/molar_mass_C*1000;
    elseif any(co2_7{x}) %ug/L
        T.CO2_converted(x)=T.CO2(x)/molar_mass_C;
    elseif any(co2_8{x}) %%sat
        T.CO2_converted(x)=T.CO2(x)*CO2_sat/100/T.Pressure(x)*T.Kh_CO2(x);
    elseif any(co2_9{x}) %mgCO2L
        T.CO2_converted(x)=T.CO2(x)*1000/molar_mass_C*(molar_mass_C/molar_mass_CO2);
    elseif any(co2_10{x}) %atm
        T.CO2_converted(x)=T.CO2(x)*T.Kh_CO2(x)*1E6;
    end
end


%% 3. Convert CH4 concentration units
%find unit type for each row in DB
units_nmol={'nmol/L','nmol L-1','nM','nmol CH4/L'}';
ch4_1=cellfun(@(c)strcmp(c,units_nmol),T.CH4_units,'UniformOutput',false);
units_uatm={'µatm','uatm','_atm','uatm at 25degC'}';
ch4_2=cellfun(@(c)strcmp(c,units_uatm),T.CH4_units,'UniformOutput',false);
units_umol={'µM','µmol/L','µmol/l','uM','_M','_mol L-1','umol/L','_mol/L','_mol CH4/L','µmol CH4/L','umolL-1'}';
ch4_3=cellfun(@(c)strcmp(c,units_umol),T.CH4_units,'UniformOutput',false);
units_mmol={'mM','mmol/L','mmol CH4/L','mmol CH4/kg'}';
ch4_4=cellfun(@(c)strcmp(c,units_mmol),T.CH4_units,'UniformOutput',false);
units_mgL={'mg/L','mg C/L','mgC/L'}';
ch4_5=cellfun(@(c)strcmp(c,units_mgL),T.CH4_units,'UniformOutput',false);
units_ugL={'ug/L','ug L-1'}';
ch4_6=cellfun(@(c)strcmp(c,units_ugL),T.CH4_units,'UniformOutput',false);
units_sat={'% sat','% saturation','%sat','%Sat'}';
ch4_7=cellfun(@(c)strcmp(c,units_sat),T.CH4_units,'UniformOutput',false);
units_mgCH4L={'mg CH4/L'}';
ch4_8=cellfun(@(c)strcmp(c,units_mgCH4L),T.CH4_units,'UniformOutput',false);
units_ugCH4L={'µg CH4/L'}';
ch4_9=cellfun(@(c)strcmp(c,units_ugCH4L),T.CH4_units,'UniformOutput',false);
units_ppb={'ppb'}';
ch4_10=cellfun(@(c)strcmp(c,units_ppb),T.CH4_units,'UniformOutput',false);

%convert all CH4 units into umol/L
T.CH4_converted=NaN(size(T,1),1);
for x=1:size(T,1)
    if any(ch4_1{x}) %nmol/L
        T.CH4_converted(x)=T.CH4(x)/1000;
    elseif any(ch4_2{x}) %uatm
        T.CH4_converted(x)=T.CH4(x)*T.Kh_CH4(x);
    elseif any(ch4_3{x}) %umol/L
        T.CH4_converted(x)=T.CH4(x);
    elseif any(ch4_4{x}) %mmol/L
        T.CH4_converted(x)=T.CH4(x)*1000;
    elseif any(ch4_5{x}) %mg/L
        T.CH4_converted(x)=T.CH4(x)/molar_mass_C*1000;
    elseif any(ch4_6{x}) %ug/L
        T.CH4_converted(x)=T.CH4(x)/molar_mass_C;
    elseif any(ch4_7{x}) %%sat
        T.CH4_converted(x)=T.CH4(x)*CH4_sat/100/T.Pressure(x)*T.Kh_CH4(x); %%%
    elseif any(ch4_8{x}) %%mgCH4L
        T.CH4_converted(x)=T.CH4(x)*1000/molar_mass_C*(molar_mass_C/molar_mass_CH4);
    elseif any(ch4_9{x}) %%ugCH4L
        T.CH4_converted(x)=T.CH4(x)/molar_mass_C*(molar_mass_C/molar_mass_CH4);
    elseif any(ch4_10{x}) %ppb
        T.CH4_converted(x)=T.CH4(x)/1000/T.Pressure(x)*T.Kh_CH4(x);
    end
end


%% 4. Convert N2O concentration units
if ismember('N2O_units',T.Properties.VariableNames)==1
    %find unit type for each row in DB
    units_nmol={'nmol/L','nmol L-1','nmol/l','nM','nmol N2O/L'}';
    n2o_1=cellfun(@(c)strcmp(c,units_nmol),T.N2O_units,'UniformOutput',false);
    units_umol={'umol L-1','µM','µmol/L','µmol/l','uM','_M','_mol L-1','umol/L','_mol/L'}';
    n2o_2=cellfun(@(c)strcmp(c,units_umol),T.N2O_units,'UniformOutput',false);
    units_ugN2OL={'µg N2O/L'}';
    n2o_3=cellfun(@(c)strcmp(c,units_ugN2OL),T.N2O_units,'UniformOutput',false);
    units_ugNL={'_g N/L'}';
    n2o_3b=cellfun(@(c)strcmp(c,units_ugNL),T.N2O_units,'UniformOutput',false);
    units_sat={'% sat','% saturation','%sat','%Sat'}';
    n2o_4=cellfun(@(c)strcmp(c,units_sat),T.N2O_units,'UniformOutput',false);

    %convert all N2O units into umol/L
    T.N2O_converted=NaN(size(T,1),1);
    for x=1:size(T,1)
        if any(n2o_1{x}) %nmol/L
            T.N2O_converted(x)=T.N2O(x)/1000;
        elseif any(n2o_2{x}) %umol
            T.N2O_converted(x)=T.N2O(x);
        elseif any(n2o_3{x}) %ugN2OL
            T.N2O_converted(x)=T.N2O(x)/molar_mass_N*(molar_mass_N/molar_mass_N2O);
        elseif any(n2o_3b{x}) %ugNL
            T.N2O_converted(x)=T.N2O(x)/molar_mass_N;
        elseif any(n2o_4{x}) %sat
            T.N2O_converted(x)=T.N2O(x)*N2O_sat/100/T.Pressure(x)*T.Kh_N2O(x);
        end
    end
end


%% 5. Convert CO2 flux units
%find unit type for each row in DB
% when not specified, it is considered that the flux is in grams of C.
units_mmolm2d={'mmol m-2 day-1','mmol/m2/d','mmol m-2 d-1','mmol C m-2 d-1','mmol m_2 d_1','mmol C/m2/d','mmol CO2 m-2 d-1','mmol CO2/m2/d','mmolm-2d-1'}';
fco2_1=cellfun(@(c)strcmp(c,units_mmolm2d),T.F_CO2_units,'UniformOutput',false);
units_mgm2h={'mg/m2/h','mg C m-2 h-1','mg m-2 h-1','mgC m-2 h-1','mgC m-2 h-1','mg m-2 hr-1','mgC/m2/h','mg m-2h-1'}';
fco2_2=cellfun(@(c)strcmp(c,units_mgm2h),T.F_CO2_units,'UniformOutput',false);
units_molm2y={'mol/m2/yr','mol/m2/y'}';
fco2_3=cellfun(@(c)strcmp(c,units_molm2y),T.F_CO2_units,'UniformOutput',false);
units_umolm2s={'umol m-2 s-1','umol/m2/s','µmol m-2 s-1','umol C/m2/s','_mol CO2 m-2 s-1','umol CO2 m-2 s-1','µmol CO2/m2/s'}';
fco2_4=cellfun(@(c)strcmp(c,units_umolm2s),T.F_CO2_units,'UniformOutput',false);
units_gm2d={'g C m-2 d-1','gC m-2 d-1','g/m2/d','g C/m2/d','g m-2 d-1','g CO2-C m-2 d-1','gC/m2d','gC/m2/d'}';
fco2_5=cellfun(@(c)strcmp(c,units_gm2d),T.F_CO2_units,'UniformOutput',false);
units_mgm2d={'mgC/m2/d','mg C m-2 d-1','mg/m2/d','mg C/m2/d','mg m-2 d-1','mgC m-2 d-1'}';
fco2_6=cellfun(@(c)strcmp(c,units_mgm2d),T.F_CO2_units,'UniformOutput',false);
units_gm2y={'g C m-2 y-1','g C/m2/y','g/m2/yr'}';
fco2_7=cellfun(@(c)strcmp(c,units_gm2y),T.F_CO2_units,'UniformOutput',false);
units_molm2h={'mol m-2 h-1'}';
fco2_8=cellfun(@(c)strcmp(c,units_molm2h),T.F_CO2_units,'UniformOutput',false);
units_gCO2m2d={'gCO2 m-2 d-1','g CO2/d','g CO2/m2/d','gCO2/m2/d'}';
fco2_9=cellfun(@(c)strcmp(c,units_gCO2m2d),T.F_CO2_units,'UniformOutput',false);
units_Mghay={'Mg C ha-1 y-1'}';
fco2_10=cellfun(@(c)strcmp(c,units_Mghay),T.F_CO2_units,'UniformOutput',false);
units_mgCO2m2d={'mg CO2/m2/d','mgCO2 m-2 d-1'}';
fco2_11=cellfun(@(c)strcmp(c,units_mgCO2m2d),T.F_CO2_units,'UniformOutput',false);
units_molm2d={'mol C/m2/d','mol/m2/d','mol CO2 m-2 d-1','mol m-2 d-1'}';
fco2_12=cellfun(@(c)strcmp(c,units_molm2d),T.F_CO2_units,'UniformOutput',false);
units_umolm2h={'umol m-2 h-1'}';
fco2_13=cellfun(@(c)strcmp(c,units_umolm2h),T.F_CO2_units,'UniformOutput',false);
units_umolm2d={'umol m-2 d-1'}';
fco2_14=cellfun(@(c)strcmp(c,units_umolm2d),T.F_CO2_units,'UniformOutput',false);
units_ugm2h={'ug C m-2 h-1','ug m-2 h-1'}';
fco2_15=cellfun(@(c)strcmp(c,units_ugm2h),T.F_CO2_units,'UniformOutput',false);
units_mmolm2h={'mmol m-2 h-1','mmol/m2/h','mmol m-2 hr-1'}';
fco2_16=cellfun(@(c)strcmp(c,units_mmolm2h),T.F_CO2_units,'UniformOutput',false);
units_umolm2min={'umol m-2 min-1'}';
fco2_17=cellfun(@(c)strcmp(c,units_umolm2min),T.F_CO2_units,'UniformOutput',false);
units_gm2h={'g/m2/h','gC m-2 h-1','g CO2-C m-2 hr-1'}';
fco2_18=cellfun(@(c)strcmp(c,units_gm2h),T.F_CO2_units,'UniformOutput',false);
units_ugm2s={'_g/m2/s','ug C m-2 s-1'}';
fco2_19=cellfun(@(c)strcmp(c,units_ugm2s),T.F_CO2_units,'UniformOutput',false);
units_nmolm2s={'nmol m-2 s-1'}';
fco2_20=cellfun(@(c)strcmp(c,units_nmolm2s),T.F_CO2_units,'UniformOutput',false);
units_gm2mon={'gC m-2 mon-1'}';
fco2_21=cellfun(@(c)strcmp(c,units_gm2mon),T.F_CO2_units,'UniformOutput',false);
units_mgCO2m2h={'mg CO2 m-2 h-1','mgCO2/m2/h'}';
fco2_22=cellfun(@(c)strcmp(c,units_mgCO2m2h),T.F_CO2_units,'UniformOutput',false);
units_nmolcm2s={'nmol cm-2 s-1'}';
fco2_23=cellfun(@(c)strcmp(c,units_nmolcm2s),T.F_CO2_units,'UniformOutput',false);
units_mmolm2s={'mmol m-2 s-1'}';
fco2_24=cellfun(@(c)strcmp(c,units_mmolm2s),T.F_CO2_units,'UniformOutput',false);

%convert all CO2 fluxes into mmol/m2/d
T.F_CO2_converted=NaN(size(T,1),1);
for x=1:size(T,1)
    if any(fco2_1{x}) %mmolm2d
        T.F_CO2_converted(x)=T.F_CO2(x);
    elseif any(fco2_2{x}) %mgm2h
        T.F_CO2_converted(x)=T.F_CO2(x)*hours_in_a_day/molar_mass_C;
    elseif any(fco2_3{x}) %molm2y
        T.F_CO2_converted(x)=T.F_CO2(x)/days_in_a_year*1000;
    elseif any(fco2_4{x}) %umolm2s
        T.F_CO2_converted(x)=T.F_CO2(x)/1000*seconds_in_a_day;
    elseif any(fco2_5{x}) %gm2d
        T.F_CO2_converted(x)=T.F_CO2(x)*1000/molar_mass_C;
    elseif any(fco2_6{x}) %mgm2d
        T.F_CO2_converted(x)=T.F_CO2(x)/molar_mass_C;
    elseif any(fco2_7{x}) %gm2y
        T.F_CO2_converted(x)=T.F_CO2(x)*1000/molar_mass_C/days_in_a_year;
    elseif any(fco2_8{x}) %molm2h
        T.F_CO2_converted(x)=T.F_CO2(x)*1000*hours_in_a_day;
    elseif any(fco2_9{x}) %gCO2m2d
        T.F_CO2_converted(x)=T.F_CO2(x)*1000/molar_mass_C*(molar_mass_C/molar_mass_CO2);
    elseif any(fco2_10{x}) %Mghay
        T.F_CO2_converted(x)=T.F_CO2(x)*milligrams_in_a_tonne/m2_in_a_ha/days_in_a_year/molar_mass_C;
    elseif any(fco2_11{x}) %mgCO2m2d
        T.F_CO2_converted(x)=T.F_CO2(x)/molar_mass_C*(molar_mass_C/molar_mass_CO2);
    elseif any(fco2_12{x}) %molm2d
        T.F_CO2_converted(x)=T.F_CO2(x)*1000;
    elseif any(fco2_13{x}) %umolm2h
        T.F_CO2_converted(x)=T.F_CO2(x)/1000*hours_in_a_day;
    elseif any(fco2_14{x}) %umolm2d
        T.F_CO2_converted(x)=T.F_CO2(x)/1000;
    elseif any(fco2_15{x}) %ugm2h
        T.F_CO2_converted(x)=T.F_CO2(x)/1000*hours_in_a_day/molar_mass_C;
    elseif any(fco2_16{x}) %mmolm2h
        T.F_CO2_converted(x)=T.F_CO2(x)*hours_in_a_day;
    elseif any(fco2_17{x}) %umolm2min
        T.F_CO2_converted(x)=T.F_CO2(x)/1000*minutes_in_a_day;
    elseif any(fco2_18{x}) %gm2h
        T.F_CO2_converted(x)=T.F_CO2(x)*1000*hours_in_a_day/molar_mass_C;
    elseif any(fco2_19{x}) %ugm2s
        T.F_CO2_converted(x)=T.F_CO2(x)/1000*seconds_in_a_day/molar_mass_C;
    elseif any(fco2_20{x}) %nmolm2s
        T.F_CO2_converted(x)=T.F_CO2(x)/1000000*seconds_in_a_day;
    elseif any(fco2_21{x}) %gm2mon
        T.F_CO2_converted(x)=T.F_CO2(x)*1000/molar_mass_C/days_in_a_month;
    elseif any(fco2_22{x}) %mgCO2m2h
        T.F_CO2_converted(x)=T.F_CO2(x)/molar_mass_C*(molar_mass_C/molar_mass_CO2)*hours_in_a_day;
    elseif any(fco2_23{x}) %nmolcm2s
        T.F_CO2_converted(x)=T.F_CO2(x)/1000000*seconds_in_a_day*10000;
    elseif any(fco2_24{x}) %mmolm2s
        T.F_CO2_converted(x)=T.F_CO2(x)*seconds_in_a_day;
    end
end


%% 6. Convert CH4 flux units
% when not specified in unit, it is considered that the flux is in grams of C.
%find unit type for each row in DB
units_mmolm2d={'mmol/m2/d','mmol m-2 d-1','mmol m_2 d_1','mmol CH4/m2/d','mmolám-2ád-1','mmolC/m2SA/d','mmol/m2/day','mmol /m2/d','mmolm-2d-1'}';
fch4_1=cellfun(@(c)strcmp(c,units_mmolm2d),T.F_CH4_units,'UniformOutput',false);
units_mgm2h={'mg/m2/h','mg C m-2 h-1','mg m-2 h-1','mgC m-2 h-1','mgC m-2 h-1','mg/(m2/h)','mgC/m2/h','mg CH4-C m-2 h-1'}';
fch4_2=cellfun(@(c)strcmp(c,units_mgm2h),T.F_CH4_units,'UniformOutput',false);
units_molm2y={'mol/m2/yr','mol/m2/y'}';
fch4_3=cellfun(@(c)strcmp(c,units_molm2y),T.F_CH4_units,'UniformOutput',false);
units_umolm2s={'umol m-2 s-1','umol/m2/s','µmol m-2 s-1','umol C/m2/s','_mol CO2 m-2 s-1','umol CO2 m-2 s-1','µmol CO2/m2/s'}';
fch4_4=cellfun(@(c)strcmp(c,units_umolm2s),T.F_CH4_units,'UniformOutput',false);
units_gm2d={'g C m-2 d-1','gC m-2 d-1','g/m2/d','g C/m2/d','g m-2 d-1','gC/m2d','gC/m2/d'}';
fch4_5=cellfun(@(c)strcmp(c,units_gm2d),T.F_CH4_units,'UniformOutput',false);
units_mgm2d={'mgC/m2/d','mg C m-2 d-1','mg/m2/d','mg C/m2/d','mg m-2 d-1','mgC m-2 d-1','mg C m_2 d_1','mg/m2/day','mg/m^2/day','mg CH4-C m-2 d-1'}';
fch4_6=cellfun(@(c)strcmp(c,units_mgm2d),T.F_CH4_units,'UniformOutput',false);
units_gm2y={'g C m-2 y-1','g C/m2/y','g/m2/yr','g C_m-2 y-1','gC/m2/y'}';
fch4_7=cellfun(@(c)strcmp(c,units_gm2y),T.F_CH4_units,'UniformOutput',false);
units_molm2h={'mol m-2 h-1'}';
fch4_8=cellfun(@(c)strcmp(c,units_molm2h),T.F_CH4_units,'UniformOutput',false);
units_mgCH4m2h={'mg CH4 /m2/h','mg CH4 * m-2 * h-1','mg CH4 m-2 h-1','mgCH4/m2/h'}';
fch4_10=cellfun(@(c)strcmp(c,units_mgCH4m2h),T.F_CH4_units,'UniformOutput',false);
units_mgCH4m2d={'mg CH4 /m2/d','mg CH4/m2/d','mgCH4/m2/d','mgCH4 m-2 d-1','mg CH4 m-2 d-1','mg CH4 m-2 day-1'};
fch4_11=cellfun(@(c)strcmp(c,units_mgCH4m2d),T.F_CH4_units,'UniformOutput',false);
units_molm2d={'mol C/m2/d','mol/m2/d','mol CH4/m2/d'}';
fch4_12=cellfun(@(c)strcmp(c,units_molm2d),T.F_CH4_units,'UniformOutput',false);
units_umolm2h={'µmol CH4/m2/h','µmol/m2/h','umol m-2 h-1'}';
fch4_13=cellfun(@(c)strcmp(c,units_umolm2h),T.F_CH4_units,'UniformOutput',false);
units_umolm2d={'umol m-2 d-1','µmol/m2/d','_mol/m2/d','_mol m_2 d_1','_mol m-2 d-1'}';
fch4_14=cellfun(@(c)strcmp(c,units_umolm2d),T.F_CH4_units,'UniformOutput',false);
units_ugm2h={'ug m-2 h-1','_g m-2 h-1','ugC m-2 h-1'}';
fch4_15=cellfun(@(c)strcmp(c,units_ugm2h),T.F_CH4_units,'UniformOutput',false);
units_mmolm2h={'mmol m-2 h-1','mmol/m2/h'}';
fch4_16=cellfun(@(c)strcmp(c,units_mmolm2h),T.F_CH4_units,'UniformOutput',false);
units_umolm2min={'umol m-2 min-1'}';
fch4_17=cellfun(@(c)strcmp(c,units_umolm2min),T.F_CH4_units,'UniformOutput',false);
units_gm2h={'g/m2/h','gC m-2 h-1'}';
fch4_18=cellfun(@(c)strcmp(c,units_gm2h),T.F_CH4_units,'UniformOutput',false);
units_mgm2min={'mg/m2/min'}';
fch4_19=cellfun(@(c)strcmp(c,units_mgm2min),T.F_CH4_units,'UniformOutput',false);
units_nmolm2s={'nmol m-2 s-1'}';
fch4_20=cellfun(@(c)strcmp(c,units_nmolm2s),T.F_CH4_units,'UniformOutput',false);
units_gm2mon={'gC m-2 mon-1'}';
fch4_21=cellfun(@(c)strcmp(c,units_gm2mon),T.F_CH4_units,'UniformOutput',false);
units_nmolcm2s={'nmol/cm2/s','nmol cm-2 s-1'}';
fch4_22=cellfun(@(c)strcmp(c,units_nmolcm2s),T.F_CH4_units,'UniformOutput',false);
units_ugCH4m2s={'ug CH4 m_2 s_1'}';
fch4_23=cellfun(@(c)strcmp(c,units_ugCH4m2s),T.F_CH4_units,'UniformOutput',false);
units_ghad={'g/ha/d','g ha-1 d-1'}';
fch4_24=cellfun(@(c)strcmp(c,units_ghad),T.F_CH4_units,'UniformOutput',false);
units_gCH4m2d={'g CH4 m-2 d-1'};
fch4_25=cellfun(@(c)strcmp(c,units_gCH4m2d),T.F_CH4_units,'UniformOutput',false);


%convert all CH4 fluxes into mmol/m2/d
% DIFFUSION
T.F_CH4_converted=NaN(size(T,1),1);
for x=1:size(T,1)
    if any(fch4_1{x}) %mmolm2d
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x);
    elseif any(fch4_2{x}) %mgm2h
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)*hours_in_a_day/molar_mass_C;
    elseif any(fch4_3{x}) %molm2y
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/days_in_a_year*1000;
    elseif any(fch4_4{x}) %umolm2s
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/1000*seconds_in_a_day;
    elseif any(fch4_5{x}) %gm2d
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)*1000/molar_mass_C;
    elseif any(fch4_6{x}) %mgm2d
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/molar_mass_C;
    elseif any(fch4_7{x}) %gm2y
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)*1000/molar_mass_C/days_in_a_year;
    elseif any(fch4_8{x}) %molm2h
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)*1000*hours_in_a_day;
    elseif any(fch4_10{x}) %mgCH4m2h
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/molar_mass_C*hours_in_a_day*(molar_mass_C/molar_mass_CH4);
    elseif any(fch4_11{x}) %mgCH4m2d
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/molar_mass_C*(molar_mass_C/molar_mass_CH4);
    elseif any(fch4_12{x}) %molm2d
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)*1000;
    elseif any(fch4_13{x}) %umolm2h
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/1000*hours_in_a_day;
    elseif any(fch4_14{x}) %umolm2d
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/1000;
    elseif any(fch4_15{x}) %ugm2h
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/1000*hours_in_a_day/molar_mass_C;
    elseif any(fch4_16{x}) %mmolm2h
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)*hours_in_a_day;
    elseif any(fch4_17{x}) %umolm2min
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/1000*minutes_in_a_day;
    elseif any(fch4_18{x}) %gm2h
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)*1000*hours_in_a_day/molar_mass_C;
    elseif any(fch4_19{x}) %mgm2min
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)*minutes_in_a_day/molar_mass_C;
    elseif any(fch4_20{x}) %nmolm2s
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/1000000*seconds_in_a_day;
    elseif any(fch4_21{x}) %gm2mon
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)*1000/molar_mass_C/days_in_a_month;
    elseif any(fch4_22{x}) %nmolcm2s
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/1000000*seconds_in_a_day*10000;
    elseif any(fch4_23{x}) %ugCH4m2s
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)/molar_mass_C*seconds_in_a_day*(molar_mass_C/molar_mass_CH4)/1000;
    elseif any(fch4_24{x}) %ghad
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)*1000/molar_mass_C/m2_in_a_ha;
    elseif any(fch4_25{x}) %gCH4m2d
        T.F_CH4_converted(x)=T.FCH4_diff_avg_val(x)*1000/molar_mass_C*(molar_mass_C/molar_mass_CH4);
    end
end

% EBULLITION
if ismember('FCH4_ebb_avg_val',T.Properties.VariableNames)==1
    T.F_CH4_ebb_converted=NaN(size(T,1),1);
    for x=1:size(T,1)
        if any(fch4_1{x}) %mmolm2d
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x);
        elseif any(fch4_2{x}) %mgm2h
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)*hours_in_a_day/molar_mass_C;
        elseif any(fch4_3{x}) %molm2y
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/days_in_a_year*1000;
        elseif any(fch4_4{x}) %umolm2s
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/1000*seconds_in_a_day;
        elseif any(fch4_5{x}) %gm2d
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)*1000/molar_mass_C;
        elseif any(fch4_6{x}) %mgm2d
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/molar_mass_C;
        elseif any(fch4_7{x}) %gm2y
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)*1000/molar_mass_C/days_in_a_year;
        elseif any(fch4_8{x}) %molm2h
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)*1000*hours_in_a_day;
        elseif any(fch4_10{x}) %mgCH4m2h
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/molar_mass_C*hours_in_a_day*(molar_mass_C/molar_mass_CH4);
        elseif any(fch4_11{x}) %mgCH4m2d
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/molar_mass_C*(molar_mass_C/molar_mass_CH4);
        elseif any(fch4_12{x}) %molm2d
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)*1000;
        elseif any(fch4_13{x}) %umolm2h
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/1000*hours_in_a_day;
        elseif any(fch4_14{x}) %umolm2d
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/1000;
        elseif any(fch4_15{x}) %ugm2h
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/1000*hours_in_a_day/molar_mass_C;
        elseif any(fch4_16{x}) %mmolm2h
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)*hours_in_a_day;
        elseif any(fch4_17{x}) %umolm2min
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/1000*minutes_in_a_day;
        elseif any(fch4_18{x}) %gm2h
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)*1000*hours_in_a_day/molar_mass_C;
        elseif any(fch4_19{x}) %mgm2min
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)*minutes_in_a_day/molar_mass_C;
        elseif any(fch4_20{x}) %nmolm2s
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/1000000*seconds_in_a_day;
        elseif any(fch4_21{x}) %gm2mon
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)*1000/molar_mass_C/days_in_a_month;
        elseif any(fch4_22{x}) %nmolcm2s
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/1000000*seconds_in_a_day*10000;
        elseif any(fch4_23{x}) %ugCH4m2s
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)/molar_mass_C*seconds_in_a_day*(molar_mass_C/molar_mass_CH4)/1000;
        elseif any(fch4_24{x}) %ghad
            T.F_CH4_ebb_converted(x)=T.FCH4_ebb_avg_val(x)*1000/molar_mass_C/m2_in_a_ha;
        end
    end
end


%% 7. Convert N2O flux units
% when not specified in unit, it is considered that the flux is in grams of N.
%find unit type for each row in DB
if ismember('FN2O_UNITS',T.Properties.VariableNames)==1
    units_umolm2d={'umol/m2/day','umol/m2/d','umol m-2 d-1','µmol/m2/d','_mol/m2/d','_mol m_2 d_1','_mol m-2 d-1','umolm-2d-1'}';
    fn2o_1=cellfun(@(c)strcmp(c,units_umolm2d),T.FN2O_UNITS,'UniformOutput',false);
    units_mgm2d={'mg/m2/day','mg N m-2 d-1','mg/m2/d','mg  m-2 d-1'}';
    fn2o_2=cellfun(@(c)strcmp(c,units_mgm2d),T.FN2O_UNITS,'UniformOutput',false);
    units_mgN2Om2d={'mg N2O m-2 day-1','mg N2O/m2/d'};
    fn2o_3=cellfun(@(c)strcmp(c,units_mgN2Om2d),T.FN2O_UNITS,'UniformOutput',false);
    units_nmolm2s={'nmol m-2 s-1'}';
    fn2o_4=cellfun(@(c)strcmp(c,units_nmolm2s),T.FN2O_UNITS,'UniformOutput',false);
    units_ugm2h={'_g N m_2 h_1','ug N m-2 h-1'}';
    fn2o_5=cellfun(@(c)strcmp(c,units_ugm2h),T.FN2O_UNITS,'UniformOutput',false);
    units_umolm2h={'umol/m2/hour','µmol/m2/h','umol m-2 h-1'}';
    fn2o_6=cellfun(@(c)strcmp(c,units_umolm2h),T.FN2O_UNITS,'UniformOutput',false);
    units_gm2y={'g/m2/yr'}';
    fn2o_7=cellfun(@(c)strcmp(c,units_gm2y),T.FN2O_UNITS,'UniformOutput',false);
    units_mmolm2d={'mmol/m2/d','mmol/m2/day'}';
    fn2o_8=cellfun(@(c)strcmp(c,units_mmolm2d),T.FN2O_UNITS,'UniformOutput',false);
    units_ugN2Om2h={'ug N2O/m2/h','µg N2O m-2 h-1'}';
    fn2o_9=cellfun(@(c)strcmp(c,units_ugN2Om2h),T.FN2O_UNITS,'UniformOutput',false);
    units_ugm2d={'ug m-2 d-1'}';
    fn2o_10=cellfun(@(c)strcmp(c,units_ugm2d),T.FN2O_UNITS,'UniformOutput',false);

    %convert all N2O fluxes into mmol/m2/d
    T.F_N2O_converted=NaN(size(T,1),1);
    for x=1:size(T,1)
        if any(fn2o_1{x}) %umolm2d
            T.F_N2O_converted(x)=T.FN2O_avg_val(x)/1000;
        elseif any(fn2o_2{x}) %mgm2d
            T.F_N2O_converted(x)=T.FN2O_avg_val(x)/molar_mass_N;
        elseif any(fn2o_3{x}) %mgN2Om2d
            T.F_N2O_converted(x)=T.FN2O_avg_val(x)/molar_mass_N*(molar_mass_N/molar_mass_N2O);
        elseif any(fn2o_4{x}) %nmolm2s
            T.F_N2O_converted(x)=T.FN2O_avg_val(x)/1000000*seconds_in_a_day;
        elseif any(fn2o_5{x}) %ugm2h
            T.F_N2O_converted(x)=T.FN2O_avg_val(x)/1000*hours_in_a_day/molar_mass_N;
        elseif any(fn2o_6{x}) %umolm2h
            T.F_N2O_converted(x)=T.FN2O_avg_val(x)/1000*hours_in_a_day;
        elseif any(fn2o_7{x}) %gm2y
            T.F_N2O_converted(x)=T.FN2O_avg_val(x)*1000/molar_mass_N/days_in_a_year;
        elseif any(fn2o_8{x}) %mmolm2d
            T.F_N2O_converted(x)=T.FN2O_avg_val(x);
        elseif any(fn2o_9{x}) %ugN2Om2h
            T.F_N2O_converted(x)=T.FN2O_avg_val(x)/1000*hours_in_a_day/molar_mass_N*(molar_mass_N/molar_mass_N2O);
        elseif any(fn2o_10{x}) %ugm2d
            T.F_N2O_converted(x)=T.FN2O_avg_val(x)/1000/molar_mass_N;
        end
    end
end


%% 8. Convert k fluxes into k600 (m/d)
if ismember('k_units',T.Properties.VariableNames)==1
    %find unit type for each row in DB
    units_md={'m/d','m d-1','meter*day^-1','m/day','m d_1'}';
    fk_1=cellfun(@(c)strcmp(c,units_md),T.k_units,'UniformOutput',false);
    units_cmh={'cm/h','cm/hr','cm h-1','cm*h-1','cm_h_1'}';
    fk_2=cellfun(@(c)strcmp(c,units_cmh),T.k_units,'UniformOutput',false);
    %Schmidt number
    T.Sc_CO2=1911-118.11*T.Temperature+3.453*T.Temperature.^2-0.0413*T.Temperature.^3;
    %conversions
    T.k600_converted=NaN(size(T,1),1);
    for x=1:size(T,1)
        if ~isnan(T.k600(x))
            T.k600_converted(x)=T.k600(x);
        elseif ~isnan(T.k(x))
            T.k600_converted(x)=T.k(x)*(600/T.Sc_CO2(x))^-0.5; %convert kCO2 to k600
        end

        if any(fk_1{x}) %m/d
            T.k600_converted(x)=T.k600_converted(x);
        elseif any(fk_2{x}) %cmh
            T.k600_converted(x)=T.k600_converted(x)/centimeters_in_a_meter*hours_in_a_day;
        end
    end
    T=movevars(T,"k600_converted",'After',"k_units");
end


%% 9. Convert Q to m3/s
if isequal(choice,'streams/rivers')
    if ismember('Discharge_mean_units',T.Properties.VariableNames)==1
        %find unit type for each row in DB
        units_m3s={'m3/s','m3 s-1','m3*s-1','m/s'}';
        fQ_1=cellfun(@(c)strcmp(c,units_m3s),T.Discharge_mean_units,'UniformOutput',false);
        units_mmy={'mm y-1'}';
        fQ_2=cellfun(@(c)strcmp(c,units_mmy),T.Discharge_mean_units,'UniformOutput',false);
        units_m3d={'m3 d-1'}';
        fQ_3=cellfun(@(c)strcmp(c,units_m3d),T.Discharge_mean_units,'UniformOutput',false);
        units_Ls={'L/s','L s-1'}';
        fQ_4=cellfun(@(c)strcmp(c,units_Ls),T.Discharge_mean_units,'UniformOutput',false);
        units_m3y={'m3 y-1','m3*yr-1'}';
        fQ_5=cellfun(@(c)strcmp(c,units_m3y),T.Discharge_mean_units,'UniformOutput',false);
        units_km3y={'km3 yr-1','km3/y'}';
        fQ_6=cellfun(@(c)strcmp(c,units_km3y),T.Discharge_mean_units,'UniformOutput',false);
        units_m3min={'m3/min'}';
        fQ_7=cellfun(@(c)strcmp(c,units_m3min),T.Discharge_mean_units,'UniformOutput',false);
        units_cfs={'cfs'}';
        fQ_8=cellfun(@(c)strcmp(c,units_cfs),T.Discharge_mean_units,'UniformOutput',false);
        units_MLd={'ML*day-1','ML/day','ML/d'}';
        fQ_9=cellfun(@(c)strcmp(c,units_MLd),T.Discharge_mean_units,'UniformOutput',false);

        %conversions
        T.Discharge_converted=NaN(size(T,1),1);
        for x=1:size(T,1)
            if any(fQ_1{x}) %m3s
                T.Discharge_converted(x)=T.Discharge_mean(x);
            elseif any(fQ_2{x}) %mmy
                T.Discharge_converted(x)=NaN;
            elseif any(fQ_3{x}) %m3d
                T.Discharge_converted(x)=T.Discharge_mean(x)/seconds_in_a_day;
            elseif any(fQ_4{x}) %Ls
                T.Discharge_converted(x)=T.Discharge_mean(x)/litres_in_a_m3;
            elseif any(fQ_5{x}) %m3y
                T.Discharge_converted(x)=T.Discharge_mean(x)/seconds_in_a_year;
            elseif any(fQ_6{x}) %km3y
                T.Discharge_converted(x)=T.Discharge_mean(x)/seconds_in_a_year*m3_in_a_km3;
            elseif any(fQ_7{x}) %m3min
                T.Discharge_converted(x)=T.Discharge_mean(x)/seconds_in_a_minute;
            elseif any(fQ_8{x}) %cfs
                T.Discharge_converted(x)=T.Discharge_mean(x)/m3_in_a_ft3;
            elseif any(fQ_9{x}) %MLd
                T.Discharge_converted(x)=T.Discharge_mean(x)/seconds_in_a_day*m3_in_a_ML;
            end
        end
    end
end


%% 10. Convert DOC concentrations in umol/L
if ismember('DOC_units',T.Properties.VariableNames)==1
    %find unit type for each row in DB
    units_umol={'_mol/L','_mol L-1','µM','_mol L_1','uM','_M','umol/L','_mol_l_1'}';
    doc_1=cellfun(@(c)strcmp(c,units_umol),T.DOC_units,'UniformOutput',false);
    units_mmol={'mM','mmol/L','mol/m3'}';
    doc_2=cellfun(@(c)strcmp(c,units_mmol),T.DOC_units,'UniformOutput',false);
    units_mgL={'mg/L','mg C/L','gC/m3','mg L-1'}';
    doc_3=cellfun(@(c)strcmp(c,units_mgL),T.DOC_units,'UniformOutput',false);

    %convert all DOC units into umol/L
    T.DOC_converted=NaN(size(T,1),1);
    for x=1:size(T,1)
        if any(doc_1{x}) %umol
            T.DOC_converted(x)=T.DOC_value(x);
        elseif any(doc_2{x}) %mmol
            T.DOC_converted(x)=T.DOC_value(x)*1000;
        elseif any(doc_3{x}) %mg/L
            T.DOC_converted(x)=T.DOC_value(x)/molar_mass_C*1000;
        end
    end
end


%% 11. Convert POC concentrations in umol/L
if ismember('POC_units',T.Properties.VariableNames)==1
    %find unit type for each row in DB
    units_umol={'_mol/L','_mol L-1','µM','_mol L_1','uM','_M','umol/L','_mol_l_1'}';
    poc_1=cellfun(@(c)strcmp(c,units_umol),T.POC_units,'UniformOutput',false);
    units_mmol={'mM','mmol/L','mol/m3'}';
    poc_2=cellfun(@(c)strcmp(c,units_mmol),T.POC_units,'UniformOutput',false);
    units_mgL={'mg/L','mg C/L','gC/m3','mg L-1'}';
    poc_3=cellfun(@(c)strcmp(c,units_mgL),T.POC_units,'UniformOutput',false);
    units_ugL={'ug/L','ug C/L'}';
    poc_4=cellfun(@(c)strcmp(c,units_ugL),T.POC_units,'UniformOutput',false);

    %convert all POC units into umol/L
    T.POC_converted=NaN(size(T,1),1);
    for x=1:size(T,1)
        if any(poc_1{x}) %umol
            T.POC_converted(x)=T.POC_value(x);
        elseif any(poc_2{x}) %mmol
            T.POC_converted(x)=T.POC_value(x)*1000;
        elseif any(poc_3{x}) %mg/L
            T.POC_converted(x)=T.POC_value(x)/molar_mass_C*1000;
        elseif any(poc_4{x}) %ug/L
            T.POC_converted(x)=T.POC_value(x)/molar_mass_C;
        end
    end

end


%% 12. Convert DIC concentrations in umol/L
if ismember('DIC_units',T.Properties.VariableNames)==1
    %find unit type for each row in DB
    units_umol={'_mol/L','_mol L-1','µM','_mol L_1','uM','_M','umol/L','µmol/L','umol/kg'}';
    dic_1=cellfun(@(c)strcmp(c,units_umol),T.DIC_units,'UniformOutput',false);
    units_mmol={'mM','mmol/L','mol/m3','mmol L-1','mM C','mmol/kg'}';
    dic_2=cellfun(@(c)strcmp(c,units_mmol),T.DIC_units,'UniformOutput',false);
    units_mgL={'mg/L','mg C/L','gC/m3'}';
    dic_3=cellfun(@(c)strcmp(c,units_mgL),T.DIC_units,'UniformOutput',false);

    %convert all DIC units into umol/L
    T.DIC_converted=NaN(size(T,1),1);
    for x=1:size(T,1)
        if any(dic_1{x}) %umol
            T.DIC_converted(x)=T.DIC_value(x);
        elseif any(dic_2{x}) %mmol
            T.DIC_converted(x)=T.DIC_value(x)*1000;
        elseif any(dic_3{x}) %mg/L
            T.DIC_converted(x)=T.DIC_value(x)/molar_mass_C*1000;
        end
    end

end


%% 13. Convert alkalinity in umol/L
if ismember('Alkalinity_units',T.Properties.VariableNames)==1
    %find unit type for each row
    units_umol={'_mol/L','_mol L-1','µM','_mol L_1','uM','_M','umol/L','µmol/L','umol/kg','µmol L-1','_eq/L','µeq/L','ueq/L'}';
    alk_1=cellfun(@(c)strcmp(c,units_umol),T.Alkalinity_units,'UniformOutput',false);
    units_mmol={'mM','mmol/L','mol/m3','mmol L-1','mM C','mmol/kg'}';
    alk_2=cellfun(@(c)strcmp(c,units_mmol),T.Alkalinity_units,'UniformOutput',false);
    units_mgL={'mg/L','mg C/L','gC/m3','mg L-1'}';
    alk_3=cellfun(@(c)strcmp(c,units_mgL),T.Alkalinity_units,'UniformOutput',false);

    %convert all DIC units into umol/L
    T.Alk_converted=NaN(size(T,1),1);
    for x=1:size(T,1)
        if any(alk_1{x}) %umol
            T.Alk_converted(x)=T.Alkalinity(x);
        elseif any(alk_2{x}) %mmol
            T.Alk_converted(x)=T.Alkalinity(x)*1000;
        elseif any(alk_3{x}) %mg/L
            T.Alk_converted(x)=T.Alkalinity(x)/molar_mass_C*1000;
        end
    end

end


%% 14. Convert EC in uS/cm
if ismember('EC_units',T.Properties.VariableNames)==1
    %find unit type for each row in DB
    units_uScm={'uS/Cm','uS/cm','µS/cm','_s/cm'}';
    ec_1=cellfun(@(c)strcmp(c,units_uScm),T.EC_units,'UniformOutput',false);
    units_mScm={'ms cm-1','mS/cm'}';
    ec_2=cellfun(@(c)strcmp(c,units_mScm),T.EC_units,'UniformOutput',false);

    %convert all EC units into uS/cm
    T.EC_converted=NaN(size(T,1),1);
    for x=1:size(T,1)
        if any(ec_1{x})
            T.EC_converted(x)=T.EC(x);
        elseif any(ec_2{x})
            T.EC_converted(x)=T.EC(x)*1000;
        end
    end

end


%% 15. Convert DO in mg/L
if ismember('DO_units',T.Properties.VariableNames)==1
    %find unit type for each row in DB
    units_percent={'% sat','%','%Sat','%sat'}';
    do_1=cellfun(@(c)strcmp(c,units_percent),T.DO_units,'UniformOutput',false);
    units_mgL={'mg/L','mg L-1','mg/l'}';
    do_2=cellfun(@(c)strcmp(c,units_mgL),T.DO_units,'UniformOutput',false);
    units_uM={'uM'}';
    do_3=cellfun(@(c)strcmp(c,units_uM),T.DO_units,'UniformOutput',false);

    %convert all DO units into mg/L
    T.DO_converted=NaN(size(T,1),1);
    for x=1:size(T,1)
        if any(do_1{x})
            O2_eq_std=exp(7.7117-1.31403*log(T.Temp_gapfill(x)+45.93)); %equilibrium O2 conc at standard pressure, mg/L
            P_watervapour=11.8571-(3840.7/T.Temp_Kelvin(x))-(216961/T.Temp_Kelvin(x)^2); %atm
            tet=0.000975-(1.426*10^(-5)*T.Temp_gapfill(x))+(6.436*10^(-8)*T.Temp_gapfill(x)^2);
            Press_atm=T.Pressure(x)/p0; %convert pressure from bars to atm
            O2_eq=O2_eq_std*Press_atm*(((1-P_watervapour/Press_atm)*(1-tet*Press_atm))/((1-P_watervapour)*(1-tet))); %equilibrium O2 conc at non-standard pressure, mg/L
            T.DO_converted(x)=O2_eq*T.DO(x)/100;
            clear O2_eq_std P_watervapour tet O2_eq
        elseif any(do_2{x})
            T.DO_converted(x)=T.DO(x);
        elseif any(do_3{x})
            T.DO_converted(x)=T.DO(x)*molar_mass_O2/1000;
            clear DOsat
        end
    end

end


Tclean=T;

disp(' ')
disp('unit conversion done')
disp(' ')

end