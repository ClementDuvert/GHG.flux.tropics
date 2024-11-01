%% Function to average out values from sites with several measurements and reduce table
% we calculate the median of measurements to avoid outlier effects.
%
function[T,Ts]=medians_per_site(T,choice)
%
% Inputs
%  - T must be a table with all variables listed as per txt file, and must have gone through 'unit_conversion' and 'climate_classification'
%  - choice is 'streams/rivers' or 'lakes/reservoirs'
%
% Outputs
%  - T is the unchanged input table
%  - Ts is the reduced table with one row per site and median values

disp(' '); disp('starting calculation of medians per site...'); disp(' ')

%find matching latitude and longitude
Loc=[T.Latitude T.Longitude];
[LatLong,ia,idx]=unique(Loc,'rows');
[counts,g]=groupcounts(Loc); %count number of repeated values
counts=counts(~isnan(g{1})); clear g %remove NaN at end

%compute dates as datenum
T.date=datenum(T.year,T.month,T.day);

%get numerical values for direct/indirect CO2 measurements, and apply one value per site as follows:
% if all measurements are direct, enter 1 (direct); if at least one measurement is indirect, enter 0 (indirect)
T.CO2_method_num=NaN(size(T,1),1); %preallocate
T.CO2_method_num(T.CO2_method_cat=='direct')=1;
T.CO2_method_num(T.CO2_method_cat=='indirect')=0;
T=movevars(T,"CO2_method_num",'After',"CO2_method_cat");
for ll=1:max(ia)
    E=T.CO2_method_num(idx==ll);
    if sum(E)<size(E,1)
        T.CO2_method_num(idx==ll)=0; %here we are conservative; if not all measurements are direct, we consider them as indirect
    elseif sum(E)==size(E,1)
        T.CO2_method_num(idx==ll)=1;
    end
    clear E
end

%calculate medians of all measurements for single sites and for a range of variables
if isequal(choice,'streams/rivers')
    R=[T.date T.Discharge_converted T.CO2_converted T.F_CO2_converted T.CH4_converted T.F_CH4_converted T.F_CH4_ebb_converted T.N2O_converted T.F_N2O_converted T.k600_converted,T.RiverWidth,T.WatershedArea,T.Temperature,T.CO2_method_num];
elseif isequal(choice,'lakes/reservoirs')
    R=[T.date T.CO2_converted T.F_CO2_converted T.CH4_converted T.F_CH4_converted T.F_CH4_ebb_converted T.N2O_converted T.F_N2O_converted T.k600_converted T.CO2_method_num T.Temperature];
end
Av=NaN(size(ia,1),size(R,2)); %preallocate
As=NaN(size(ia,1),size(R,2));
for d=1:size(R,2)
    Av(:,d)=accumarray(idx,R(:,d),[],@(x)median(x,'omitnan')); %calculate medians while removing NaNs
    As(:,d)=accumarray(idx,R(:,d),[],@(x)sum(~isnan(x))); %calculate number of elements for each median calculation
end

%create a smaller table with averaged variables
if isequal(choice,'streams/rivers')
    Nb_CO2=As(:,3);
    Nb_CH4=As(:,5);
    Nb_N2O=As(:,8);
    Nb_FCO2=As(:,4);
    Nb_FCH4diff=As(:,6);
    Nb_FCH4ebb=As(:,7);
    Nb_FN2O=As(:,9);
    Mid_date=Av(:,1);
    Discharge_converted=Av(:,2);
    CO2_converted=Av(:,3);
    F_CO2_converted=Av(:,4);
    CH4_converted=Av(:,5);
    F_CH4_diff_converted=Av(:,6);
    F_CH4_ebb_converted=Av(:,7);
    N2O_converted=Av(:,8);
    F_N2O_converted=Av(:,9);
    k600_converted=Av(:,10);
    Width=Av(:,11);
    Catchment_area=Av(:,12);
    Temperature=Av(:,13);
    CO2_direct_indirect=Av(:,14);
    Ts=table(Mid_date,LatLong,Width,Catchment_area,Discharge_converted,CO2_converted,Nb_CO2,F_CO2_converted,Nb_FCO2,CH4_converted,Nb_CH4,F_CH4_diff_converted,Nb_FCH4diff,F_CH4_ebb_converted,Nb_FCH4ebb,N2O_converted,Nb_N2O,F_N2O_converted,Nb_FN2O,k600_converted,Temperature,CO2_direct_indirect);
elseif isequal(choice,'lakes/reservoirs')
    Nb_CO2=As(:,2);
    Nb_CH4=As(:,4);
    Nb_N2O=As(:,7);
    Nb_FCO2=As(:,3);
    Nb_FCH4diff=As(:,5);
    Nb_FCH4ebb=As(:,6);
    Nb_FN2O=As(:,8);
    Mid_date=Av(:,1);
    CO2_converted=Av(:,2);
    F_CO2_converted=Av(:,3);
    CH4_converted=Av(:,4);
    F_CH4_diff_converted=Av(:,5);
    F_CH4_ebb_converted=Av(:,6);
    N2O_converted=Av(:,7);
    F_N2O_converted=Av(:,8);
    k600_converted=Av(:,9);
    CO2_direct_indirect=Av(:,10);
    Temperature=Av(:,11);
    Ts=table(Mid_date,LatLong,CO2_converted,Nb_CO2,F_CO2_converted,Nb_FCO2,CH4_converted,Nb_CH4,F_CH4_diff_converted,Nb_FCH4diff,F_CH4_ebb_converted,Nb_FCH4ebb,N2O_converted,Nb_N2O,F_N2O_converted,Nb_FN2O,k600_converted,CO2_direct_indirect,Temperature);
end

%add other variables to new table; requires reducing/reordering original table
S=T(ia,:);
Ts.Ref=S.Ref; Ts=movevars(Ts,"Ref",'Before',"Mid_date");
Ts.SiteName=S.SiteName; Ts=movevars(Ts,"SiteName",'Before',"Mid_date");
Ts.Climate=S.Climate; Ts=movevars(Ts,"Climate",'After',"LatLong");
Ts.KG=S.KG; Ts=movevars(Ts,"KG",'After',"Climate");
Ts.Elevation=S.Elevation; Ts=movevars(Ts,"Elevation",'After',"Climate");
if isequal(choice,'streams/rivers')
    Ts.Slope=S.Slope; Ts=movevars(Ts,"Slope",'After',"Mid_date");
    Ts.StrahlerOrder=S.StrahlerOrder; Ts=movevars(Ts,"StrahlerOrder",'Before',"CO2_converted");
elseif isequal(choice,'lakes/reservoirs')
    Ts.Area=S.SurfaceArea_m2./1E6; Ts=movevars(Ts,"Area",'After',"KG"); %km2 (surface area is in m2)
    Ts.Volume=S.Volume; Ts=movevars(Ts,"Volume",'After',"Area");
    Ts.MeanDepth=S.MeanDepth; Ts=movevars(Ts,"MeanDepth",'After',"Volume");
    Ts.MaxDepth=S.MaxDepth; Ts=movevars(Ts,"MaxDepth",'After',"MeanDepth");
    Ts.LakeType=S.LakeType;
end

Ts=Ts(~isnan(Ts.LatLong(:,1)),:); %remove NaNs at end

% trim off dataset to keep only (sub)tropics
Ts=Ts(Ts.LatLong(:,1)<34.00001 & Ts.LatLong(:,1)>-34.00001,:);

disp(' ')
disp('calculation of medians per site done')
disp(' ')

end


