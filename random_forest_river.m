%% Function to create a RF model stream/river GHG concentration data
%
function[allModels,best_RF_models,allFeatureImportances,best_R2,allR2]=random_forest_river(T)
%
% Input: T is a reduced table that has been processed using 'unit_conversion' and 'averaging_per_site'
% needs 'swtest.m'
%
% Output:
% 'allModels' is a m-by-n cell (m is the nb of runs; n is the nb of objective evaluations) that contains all models
% 'allFeatureImportances' is a m-by-n cell with feature importances for each predictor based on each model
% 'best_RF_models' is a 1x3 cell that contains the best three models for each gas


%% 1. Prepare and transform data as required
disp(' '); disp('starting random forest modelling...'); disp(' ')

normality_threshold=0.05; %threshold for Shapiro-Wilk p-value
numRuns=20;
nb_obj_eval=25;

% We use the concentration data as response variables
gas={'CO2_converted','CH4_converted','N2O_converted'};

%preallocate
allR2=zeros(numRuns,length(gas));
best_R2=-Inf(1,length(gas));
allFeatureImportances=cell(numRuns,length(gas));

for g=1:length(gas)
    disp(' '); disp(['starting RF model for ',gas{g}]); disp(' ')

    %this is for C and real flux
    all_variables={gas{g},'rainfall_BasAt','airTemperature_BasAt','seasonality','slope_consolidated', ...
        'humanFootprint_BasAt','soilOrganicCarbon_BasAt','discharge','groundwaterDepth_BasAt', ...
        'netPrimaryProductivity_Modis','peatlandCover_Modis'};

    predictors={'rainfall_BasAt', 'airTemperature_BasAt','seasonality','slope_consolidated', ...
        'humanFootprint_BasAt','soilOrganicCarbon_BasAt','discharge','groundwaterDepth_BasAt',...
        'netPrimaryProductivity_Modis','peatlandCover_Modis'};

    new_names={'rainfall amount','temperature','rainfall seasonality','channel slope',...
        'human footprint','soil organic carbon',...
        'discharge','groundwater depth','net primary productivity','peatland extent'};

    % response variable: log-transform if not normally distributed
    R_log=T(~isnan(T.(gas{g})),:);
    p_data=R_log.(gas{g});
    minValue=min(p_data); c=max(1,abs(minValue)+1);
    %perform Shapiro-Wilk test for normality
    [~,p,~]=swtest(p_data,normality_threshold);
    if p<normality_threshold
        %apply log transformation
        disp(['log-transforming variable: ' gas{g}])
        if any(p_data<=0)
            disp(['we found zero or negative values in the response variable ' gas{g}])
            R_log.(gas{g})=log10(R_log.(gas{g})+c);
        else
            R_log.(gas{g})=log10(R_log.((gas{g})));
        end
    end

    % predictor variables: log-transform the 3 highly skewed variables
    R_log_pred=R_log;
    if ismember('slope_consolidated',new_names)
        if any(R_log.slope_consolidated<=0)
            R_log_pred.slope_consolidated=log10(R_log.slope_consolidated+abs(min(R_log.slope_consolidated(R_log.slope_consolidated<=0)))+1); % shift to avoid log of non-positive values
        else
            R_log_pred.slope_consolidated=log10(R_log.slope_consolidated);
        end
    end
    if any(R_log.discharge<=0)
        R_log_pred.discharge=log10(R_log.discharge+abs(min(R_log.discharge(R_log.discharge<=0)))+1); % shift to avoid log of non-positive values
    else
        R_log_pred.discharge=log10(R_log.discharge);
    end
    if any(R_log.groundwaterDepth_BasAt<=0)
        R_log_pred.groundwaterDepth_BasAt=log10(R_log.groundwaterDepth_BasAt+abs(min(R_log.groundwaterDepth_BasAt(R_log.groundwaterDepth_BasAt<=0)))+1); % shift to avoid log of non-positive values
    else
        R_log_pred.groundwaterDepth_BasAt=log10(R_log.groundwaterDepth_BasAt);
    end

    % predictor variables: scale all variables
    R_norm=R_log_pred;
    for i=1:length(predictors)
        p_data=R_norm.(predictors{i});
        p_data_cleaned=p_data(~isnan(p_data));
        p_data_cleaned(p_data_cleaned==-999)=5;
        mean_val=mean(p_data_cleaned);
        std_val=std(p_data_cleaned);
        if std_val==0  %avoid division by zero
            std_val=1;
        end
        scaled_data=(p_data_cleaned-mean_val)/std_val;
        R_norm.(predictors{i})(~isnan(p_data))=scaled_data;
        clear mean_val std_val p_data_cleaned p_data scaled_data
    end

    % keep only variables of interest and remove NaNs
    R2=R_norm(:,all_variables);
    R_norm_cleaned=R2;
    for i=1:length(all_variables)
        R_norm_cleaned=R_norm_cleaned(~isnan(R_norm_cleaned.(all_variables{i})),:);
    end

    % separate predictors and response
    R_norm_predictors=R_norm_cleaned(:,predictors);
    R_norm_predictors.Properties.VariableNames=new_names;

    %correlation plot
    numX=numel(predictors);
    corrMatrix=corr(table2array(R_norm_predictors),'Rows','complete');
    nColors=100;
    cmap=brewermap(nColors,'BrBG');
    getColorIndex=@(val) round((val+1)/2*(nColors-1))+1;
    getColor=@(val) cmap(getColorIndex(val),:);
    figure('Position',[100,100,1000,400]);
    for i=1:numX
        for j=1:i-1 % lower triangle only, excluding diagonal
            value=corrMatrix(i,j);
            color=getColor(value);
            rectangle('Position',[j-0.5,numX-i+0.5,1,1], ...
                'FaceColor',color,'EdgeColor','none'); %if we want to add some alpha: 'FaceColor',[color,alpha]
        end
    end
    set(gca,'XTick',1:numX-1,'XTickLabel',new_names(1:numX-1), ...
        'YTick',1:numX-1,'YTickLabel',flip(new_names(2:numX)), ...
        'XTickLabelRotation',90,'FontSize',10,'Color','w');
    axis square;
    %add correlation values as text
    [x,y]=meshgrid(1:numX, 1:numX);
    for i = 1:numX
        for j = 1:i-1 % Lower triangle only
            text(x(i,j),numX-y(i,j)+1,sprintf('%.2f',corrMatrix(i,j)), ...
                'HorizontalAlignment','center','FontSize',8,'Color','k');
        end
    end
    box on
    pause(0.01)


    %% 2. Run RF model
    %parameterisation of RF
    hyperparams=[
        optimizableVariable('NumLearningCycles',[300,1600],'Type','integer'), ...
        optimizableVariable('MinLeafSize',[2,25],'Type','integer'), ...
        optimizableVariable('MaxNumSplits',[10,100],'Type','integer'), ...
        optimizableVariable('NumVariablesToSample',[1,length(predictors)],'Type','integer')
        ];

    %define optimisation options
    hyperoptoptions=struct('AcquisitionFunctionName','expected-improvement-plus',...
        'KFold',5,'ShowPlots',false,'MaxObjectiveEvaluations',nb_obj_eval);  %5-fold cross-validation

    %train the model with Bayesian optimisation
    for i=1:numRuns
        rng('shuffle'); %to use a different seed each time

        %model
        results=fitrensemble(R_norm_cleaned,gas{g},'Method','Bag',...
            'OptimizeHyperparameters',hyperparams,...
            'HyperparameterOptimizationOptions',hyperoptoptions);

        %store model and feature importance values
        RF_model=results.Trained{1};
        allModels{i,g}=RF_model;
        allFeatureImportances{i,g}=predictorImportance(RF_model);

        % calculate R-squared
        y_true=R_norm_cleaned.(gas{g});
        y_pred=predict(RF_model,R_norm_cleaned);
        SS_tot=sum((y_true-mean(y_true)).^2);
        SS_res=sum((y_true-y_pred).^2);
        allR2(i,g)=1-(SS_res/SS_tot);

        % check if current model is the best one
        if allR2(i,g)>best_R2(g)
            best_R2(g)=allR2(i,g);
            best_RF_models{g}=RF_model; %store best model
        end
    end

end

disp(' ')
disp('random forest modelling done')
disp(' ')


end
