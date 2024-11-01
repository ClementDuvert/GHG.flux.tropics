%% Script to upload GHG databases and call different functions (convert units, generate medians, calculate fluxes etc.)
% requires at least one txt file in the working directory:
% 'Dataset_rivers.txt' and/or 'Dataset_lakes_ponds_reservoirs.txt'

% initialise
clearvars
%addpath('/Users/ClementDuvert/...') %add directory where all functions are stored
%addpath('/Users/ClementDuvert/...') %add directory where datasets are stored
%cd('/Users/ClementDuvert/...') %add working directory

% upload data for a chosen ecosystem type
c=questdlg('Which ecosystem type?','Ecosystem type','streams/rivers','lakes/reservoirs','streams/rivers');
switch c
    case 'streams/rivers'
        T=readtable('Dataset_rivers.txt','ReadVariableNames',true,'NumHeaderLines',0,'Delimiter','\t','TreatAsEmpty',["NA","NR"]);
    case 'lakes/reservoirs'
        T=readtable('Dataset_lakes_ponds_reservoirs.txt','ReadVariableNames',true,'NumHeaderLines',0,'Delimiter','\t','TreatAsEmpty',["NA","NR"]);
end

% convert units
[T_conv]=unit_conversion(T,c);

% add a column with 5 climate classes
[T_conv_climate]=climate_classification(T_conv);

% calculate medians for sites with multiple GHG flux/concentration values
[T_conv_climate,T_reduced]=medians_per_site(T_conv_climate,c);

% calculate weighted means for each lake, pond, and reservoir (requires the correct HydroLAKES input files, see documentation in function)
if isequal(c,'lakes/reservoirs')
    T_per_system=weighted_means_per_lake;
end

% calculate fluxes using bootstrapping approach + upscaling
if isequal(c,'streams/rivers')
    [T_reduced,arealFlux_CO2,scaledFlux_CO2,scaledFlux_eq_CO2,arealFlux_CH4,...
    scaledFlux_CH4,scaledFlux_eq_CH4,arealFlux_N2O,scaledFlux_N2O,...
    scaledFlux_eq_N2O]=flux_calculation_rivers(T_reduced);
elseif isequal(c,'lakes/reservoirs')
    [scaledFlux_eq_CO2_lake,scaledFlux_eq_CO2_reservoir,...
    scaledFlux_eq_CH4_lake,scaledFlux_eq_CH4_reservoir,...
    scaledFlux_eq_N2O_lake,scaledFlux_eq_N2O_reservoir]=flux_calculation_lakes(T_per_system);
end

% run random forest model for river GHG concentration data
if isequal(c,'streams/rivers')
    [allModels,best_RF_models,allFeatureImportances]=random_forest_river(T_reduced);
end


