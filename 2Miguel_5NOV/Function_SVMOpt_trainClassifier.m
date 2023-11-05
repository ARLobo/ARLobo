function [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)

% % % Importar para o Workspace os vectores das fiadas BS a treinar
PS_BS_Treinar = readtable("E:\BkSctt_Research\1_MatLab_MLClassif\1_ML_BkSctt Tests\Dados2Treino\EM2040_SH\TrainData_LinMedian_EM2040_Single_BSModel.txt");
% % SB_APL_Treinar = readtable("E:\BkSctt_Research\1_MatLab_MLClassif\2_DadosAPL\ASCII_Files\DadosDeTreino\SB_APLData_2TrainClassf.txt");

% Extract predictors and response
% This code processes the data into the right shape for training the model.
inputTable = PS_BS_Treinar;
predictorNames = {'ANG_MEDIO', 'BS_MEDIO', 'INTERC', 'SLOPE'};
predictors = inputTable(:, predictorNames);
response = inputTable.CLASSE;
isCategoricalPredictor = [false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateSVM(...
    'KernelFunction', 'linear', ...
    'PolynomialOrder', [], ...
    'KernelScale', 1, ...
    'BoxConstraint', 1.73380193696458, ...
    'Standardize', true);
classificationSVM = fitcecoc(...
    predictors, ...
    response, ...
    'Learners', template, ...
    'Coding', 'onevsone', ...
    'ClassNames', categorical({'AF'; 'AG'; 'AM'; 'AMF'; 'AMG'}));

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
svmPredictFcn = @(x) predict(classificationSVM, x);
trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.RequiredVariables = {'ANG_MEDIO', 'BS_MEDIO', 'INTERC', 'SLOPE'};
trainedClassifier.ClassificationSVM = classificationSVM;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2020b.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = PS_BS_Treinar;
predictorNames = {'ANG_MEDIO', 'BS_MEDIO', 'INTERC', 'SLOPE'};
predictors = inputTable(:, predictorNames);
response = inputTable.CLASSE;
isCategoricalPredictor = [false, false, false, false];

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationSVM, 'KFold', 5);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
