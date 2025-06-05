%% Function to assign catchment properties to each site in raw database
%
function[T]=catchment_properties_rivers(T)
%
% Inputs
% T is a table with raw river data
% The following files are needed in the input dataset folder:
% 'river_properties_BasinATLAS_2025.csv
% 'river_properties_GRADES_2025.csv
% 'river_properties_H90_2025.csv
% 'river_properties_npp_peat_2025.csv'
%
% Outputs
% new table T with same number of rows and an increased number of columns.

disp(' '); disp('starting addition of catchment properties...'); disp(' ')


%% 1. Import landscape properties from BasinATLAS, MODIS, Hydrography90m and GRADES
BasAt=readtable('river_properties_BasinATLAS_2025.csv');
BasAt=sortrows(BasAt,'site_numeric'); %sort rows according to site indices
Modis=readtable('river_properties_npp_peat_2025.csv');
Modis=sortrows(Modis,'site_numeric'); %sort rows according to site indices
for i=1:size(Modis,1) %add lat long to Modis
    idx=BasAt.site_numeric==Modis.site_numeric(i);
    Modis.latitude(i)=BasAt.latitude_reported(idx);
    Modis.longitude(i)=BasAt.longitude_reported(idx);
    clear idx
end
disp('BasinATLAS data imported')

H90=readtable('river_properties_H90_2025.csv');
H90=sortrows(H90,'site_numeric'); %sort rows according to site indices
GRADES=readtable('river_properties_GRADES_2025.csv');
GRADES=sortrows(GRADES,'site_numeric'); %sort rows according to site indices
for i=1:size(GRADES,1) %add lat long to GRADES
    idx=BasAt.site_numeric==GRADES.site_numeric(i);
    GRADES.latitude(i)=BasAt.latitude_reported(idx);
    GRADES.longitude(i)=BasAt.longitude_reported(idx);
    clear idx
end

disp('GRADES and H90 data imported')

%% 2. Add properties to each row in the dataset
tol=1e-5;

T.humanFootprint_BasAt=NaN(size(T,1),1);
T.soilOrganicCarbon_BasAt=NaN(size(T,1),1);
T.airTemperature_BasAt=NaN(size(T,1),1);
T.rainfall_BasAt=NaN(size(T,1),1);
T.runoff_BasAt=NaN(size(T,1),1);
T.wetlandExtent_BasAt=NaN(size(T,1),1);
T.forestExtent_BasAt=NaN(size(T,1),1);
T.aridity_BasAt=NaN(size(T,1),1);
T.groundwaterDepth_BasAt=NaN(size(T,1),1);
T.P1=NaN(size(T,1),1); T.P2=NaN(size(T,1),1); T.P3=NaN(size(T,1),1);
T.P4=NaN(size(T,1),1); T.P5=NaN(size(T,1),1); T.P6=NaN(size(T,1),1);
T.P7=NaN(size(T,1),1); T.P8=NaN(size(T,1),1); T.P9=NaN(size(T,1),1);
T.P10=NaN(size(T,1),1); T.P11=NaN(size(T,1),1); T.P12=NaN(size(T,1),1);
T.netPrimaryProductivity_Modis=NaN(size(T,1),1);
T.peatlandCover_Modis=NaN(size(T,1),1);
T.slope_H90=NaN(size(T,1),1);
T.elevation_H90=NaN(size(T,1),1);
T.streamOrder_H90=NaN(size(T,1),1);
T.catchmentArea_H90=NaN(size(T,1),1);
T.slope_GRADES=NaN(size(T,1),1);
T.catchmentArea_GRADES=NaN(size(T,1),1);
T.streamOrder_GRADES=NaN(size(T,1),1);

for e=1:size(T,1)
    idx_BasAt=abs(BasAt.longitude_reported-T.Longitude(e))<tol & abs(BasAt.latitude_reported-T.Latitude(e))<tol;
    %idx_BasAt=BasAt.longitude_reported==T.Longitude(e) & BasAt.latitude_reported==T.Latitude(e);
    idx_ones=find(idx_BasAt==1);
    if numel(idx_ones)>1
        idx_BasAt(idx_ones(2:end))=0;
    end

    if any(idx_BasAt)
        T.humanFootprint_BasAt(e)=BasAt.human_footprint_up_09(idx_BasAt);
        T.soilOrganicCarbon_BasAt(e)=BasAt.soil_org_carbon_up_tonneshectare(idx_BasAt);
        T.airTemperature_BasAt(e)=BasAt.air_temp_up_avg_celsius(idx_BasAt)/10;
        T.rainfall_BasAt(e)=BasAt.precip_up_mm(idx_BasAt);
        T.runoff_BasAt(e)=BasAt.runoff_mm(idx_BasAt);
        T.wetlandExtent_BasAt(e)=BasAt.wetland_grouped2_up_per(idx_BasAt);
        T.forestExtent_BasAt(e)=BasAt.forest_cover_up_per(idx_BasAt);
        T.aridity_BasAt(e)=BasAt.aridity_index_up(idx_BasAt);
        T.groundwaterDepth_BasAt(e)=BasAt.gw_table_cm(idx_BasAt);
        T.P1(e)=BasAt.precip_sub_s01(idx_BasAt); T.P2(e)=BasAt.precip_sub_s02(idx_BasAt); T.P3(e)=BasAt.precip_sub_s03(idx_BasAt);
        T.P4(e)=BasAt.precip_sub_s04(idx_BasAt); T.P5(e)=BasAt.precip_sub_s05(idx_BasAt); T.P6(e)=BasAt.precip_sub_s06(idx_BasAt);
        T.P7(e)=BasAt.precip_sub_s07(idx_BasAt); T.P8(e)=BasAt.precip_sub_s08(idx_BasAt); T.P9(e)=BasAt.precip_sub_s09(idx_BasAt);
        T.P10(e)=BasAt.precip_sub_s10(idx_BasAt); T.P11(e)=BasAt.precip_sub_s11(idx_BasAt); T.P12(e)=BasAt.precip_sub_s12(idx_BasAt);
    end

    if ~any(idx_BasAt)
    warning('No match for site %d (%f, %f)',e,T.Latitude(e),T.Longitude(e));
elseif sum(idx_BasAt)>1
    warning('Multiple matches (%d) for site %d (%f, %f)',sum(idx_BasAt),e,T.Latitude(e),T.Longitude(e));
    end


    %idx_Modis=Modis.longitude==T.Longitude(e) & Modis.latitude==T.Latitude(e);
    idx_Modis=abs(Modis.longitude-T.Longitude(e))<tol & abs(Modis.latitude-T.Latitude(e))<tol;
    idx_ones=find(idx_Modis==1);
    if numel(idx_ones)>1
        idx_Modis(idx_ones(2:end))=0;
    end

    if any(idx_Modis)
        T.netPrimaryProductivity_Modis(e)=Modis.NPP_yr(idx_Modis);
        T.peatlandCover_Modis(e)=Modis.peatland_cover(idx_Modis);
    end

    %idx_H90=H90.longitude_reported==T.Longitude(e) & H90.latitude_reported==T.Latitude(e);
    idx_H90=abs(H90.longitude_reported-T.Longitude(e))<tol & abs(H90.latitude_reported-T.Latitude(e))<tol;
    idx_ones=find(idx_H90==1);
    if numel(idx_ones)>1
        idx_H90(idx_ones(2:end))=0;
    end

    if any (idx_H90)
        T.slope_H90(e)=H90.slope_channel(idx_H90);
        T.elevation_H90(e)=H90.elevation_m(idx_H90);
        T.streamOrder_H90(e)=H90.order_snapped_coords(idx_H90);
        T.catchmentArea_H90(e)=H90.catchment_area_snapped(idx_H90); %km2
    end

    idx_GRADES=abs(GRADES.longitude-T.Longitude(e))<tol & abs(GRADES.latitude-T.Latitude(e))<tol;
    %idx_GRADES=GRADES.longitude==T.Longitude(e) & GRADES.latitude==T.Latitude(e);
    idx_ones=find(idx_GRADES==1);
    if numel(idx_ones)>1
        idx_GRADES(idx_ones(2:end))=0;
    end

    if any(idx_GRADES)
        T.slope_GRADES(e)=GRADES.slope_grades(idx_GRADES);
        T.catchmentArea_GRADES(e)=GRADES.upstream_area_grades(idx_GRADES)/1E4;
        T.streamOrder_GRADES(e)=GRADES.order_grades(idx_GRADES);
    end

    clear idx_BasAt idx_Modis idx_H90 idx_GRADES
end

disp('all properties added')

%further adjustments and consolidation
T.shifted_streamOrder_H90=T.streamOrder_H90-1;
T.shifted_streamOrder_H90(T.shifted_streamOrder_H90==0)=1;
T.streamorder_consolidated=T.shifted_streamOrder_H90;
%correct stream orders with actual order if reported
T.streamorder_consolidated(~isnan(T.StrahlerOrder))=T.StrahlerOrder(~isnan(T.StrahlerOrder)); 
%worst case, we use the GRADES stream order
T.streamorder_consolidated(isnan(T.streamorder_consolidated))=T.streamOrder_GRADES(isnan(T.streamorder_consolidated)); 

disp('stream order consolidated')

%choice of slope
slope=questdlg('which slope product?','SLOPE','GRADES','H90','GRADES');

switch slope
    case 'GRADES'
        T.slope_GRADES(T.slope_GRADES==0)=5E-5; %assign very low values to sites with slope 0
        T.slope_consolidated=T.slope_GRADES;
        T.slope_consolidated(~isnan(T.Slope))=T.Slope(~isnan(T.Slope)); %prefer actual slope value if known
    case 'H90'
        T.slope_consolidated=T.slope_H90; %we prefer H90 slope
        idx=T.slope_consolidated==0 | isnan(T.slope_consolidated);
        T.slope_consolidated(idx)=T.slope_GRADES(idx); %if H90 is 0 or NaN, we choose GRADES slope
        T.slope_consolidated(T.slope_consolidated==0)=5E-5; %if GRADES slope is 0 we assign a low value
        T.slope_consolidated(~isnan(T.Slope))=T.Slope(~isnan(T.Slope)); %prefer actual slope value if known
end

disp('slope consolidated')

% calculate discharge and flow velocity based on BasinATLAS (runoff) and H90 (area)
T.discharge=NaN(size(T,1),1); T.velocity=NaN(size(T,1),1);
valid_idx=~isnan(T.runoff_BasAt) & ~isnan(T.catchmentArea_H90);
T.discharge(valid_idx)=(T.runoff_BasAt(valid_idx)/1E3).*(T.catchmentArea_H90(valid_idx)*1E6)/(3600*24*365); %discharge in m3/s
T.velocity(valid_idx)=exp(0.12*log(T.discharge(valid_idx))-1.06); %m/s  %based on relationship in Liu et al. 2022 PNAS ('log' is natural logarithm)

disp('discharge and velocity calculated')

%calculate rainfall seasonality as per Feng et al. 2013
T.Seasonality=NaN(height(T),1);
pmax=max(T.rainfall_BasAt);
T.Pyear=T.P1+T.P2+T.P3+T.P4+T.P5+T.P6+T.P7+T.P8+T.P9+T.P10+T.P11+T.P12; %we don't use T.rainfall_BasAt because there can be a small discrepancy between the two
for e=1:size(T,1)
    pm=NaN(1,12); D_month=zeros(1,12);
    for month=1:12
        tx=sprintf('P%.0f',month);
        if T.Pyear(e)>0
            pm(month)=T.(tx)(e)/T.Pyear(e); %monthly probability distribution
        else
            pm(month)=NaN;
        end

        if pm(month)>0
            D_month(month)=pm(month)*log2(pm(month)/(1/12)); %relative entropy
        else
            D_month(month)=0;
        end
    end
    T.Seasonality(e)=sum(D_month).*T.Pyear(e)./pmax; %rainfall seasonality
end


disp('seasonality index added')

disp(' ')
disp('catchment properties added')
disp(' ')

