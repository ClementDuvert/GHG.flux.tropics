%% Script to upload GHG databases and call different functions (convert units, generate medians, calculate fluxes etc.)
%
% requires at least one txt file in the working directory:
% 'dataset_rivers_raw.txt' and/or 'dataset_lakes_raw.txt'
%
% initialise
clearvars
%addpath('/Users/ClementDuvert/...') %add directory where functions are stored
%addpath('/Users/ClementDuvert/...') %add directory where datasets are stored
%cd('/Users/ClementDuvert/...') %add working directory

% upload data for a chosen ecosystem type
c=questdlg('Which ecosystem type?','Ecosystem type','streams/rivers','lakes/reservoirs','streams/rivers');
switch c
    case 'streams/rivers'
        T=readtable('dataset_rivers_raw.txt','ReadVariableNames',true,'NumHeaderLines',0,'Delimiter','\t','TreatAsEmpty',["NA","NR"]);
    case 'lakes/reservoirs'
        T=readtable('dataset_lakes_raw.txt','ReadVariableNames',true,'NumHeaderLines',0,'Delimiter','\t','TreatAsEmpty',["NA","NR"]);
end

% add a column with 5 climate classes
[T_climate]=climate_classification(T);

% add columns with catchment properties
if isequal(c,'streams/rivers')
    [T_properties]=catchment_properties_rivers(T_climate);
end

% convert units
if isequal(c,'streams/rivers')
    [T_conv]=unit_conversion(T_properties,c);
elseif isequal(c,'lakes/reservoirs')
    [T_conv]=unit_conversion(T_climate,c);
end

% calculate means/medians for sites/systems with multiple GHG flux/concentration values - merge sites that are <5km away from each other in large rivers
[T_conv,T_reduced]=averaging_per_site(T_conv,c);

% calculate fluxes using bootstrapping approach + upscaling
if isequal(c,'streams/rivers')
    [T_reduced,errorBar_CO2,errorBar_CH4,errorBar_N2O,arealFlux_CO2,scaledFlux_CO2,scaledFlux_eq_CO2,arealFlux_CH4,...
    scaledFlux_CH4,scaledFlux_eq_CH4,arealFlux_N2O,scaledFlux_N2O,...
    scaledFlux_eq_N2O]=flux_calculation_rivers(T_reduced);
elseif isequal(c,'lakes/reservoirs')
    [T_reduced,scaledFlux_eq_CO2_lake,scaledFlux_eq_CO2_reservoir,...
    scaledFlux_eq_CH4_lake,scaledFlux_eq_CH4_reservoir,...
    scaledFlux_eq_N2O_lake,scaledFlux_eq_N2O_reservoir]=flux_calculation_lakes(T_reduced);
end

% run random forest model for river GHG concentration data
if isequal(c,'streams/rivers')
    [allModels,best_RF_models,allFeatureImportances,best_R2,allR2]=random_forest_river(T_reduced);
end


