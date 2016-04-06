function [ x ] = WienerGapFixer(x_gap, i_gapStart, gapLength, traindelta)
%WINERGAPFIXER Summary of this function goes here
%   Detailed explanation goes here

trainL = ceil(traindelta/2); 

i_trainLeft = (i_gapStart>traindelta)*(i_gapStart-traindelta) + (i_gapStart<=traindelta);
i_trainRigth = (length(x_gap)- (i_gapStart+gapLength) > traindelta)*(i_gapStart+gapLength+traindelta) + (length(x_gap)- (i_gapStart+gapLength) <= traindelta)*length(x_gap);


data = cell(1,2); 
data{1} = x_gap(i_trainLeft:i_gapStart-1);
data{2} = flipud(x_gap(i_gapStart+gapLength:i_trainRigth));

trainL = min(length(data{1}),trainL);
trainL = min(length(data{2}),trainL);

constructed = zeros(gapLength,2);
way = [1 -1];
for n = 1:2
    
    trainData = data{n}; 
    mu = mean(trainData);
    trainData = trainData - mu; 
    
    [a,g] = lpc(trainData, trainL);

    a(find(isnan(a))) = 0;
    
    predict = filter([0 -a(2:end)],1,[trainData(1:end-1); 0]);
    errorCor = max(predict(end)/trainData(end),1.022);
    
    for index = 1:gapLength;
        predict = filter([0 -a(2:end)],1,[trainData; 0]);
        trainData = [trainData; predict(end)*errorCor];
        constructed(way(n)*index+(((way(n)-1)/2)*(-end-1)), n) = predict(end)*errorCor;        
    end
    
    constructed(:, n) = constructed(:, n) + mu;
end

prop = linspace(1,0, gapLength)' ;

constructed(isnan(constructed)) = 0;

G = [ prop  flipud(prop)];
V = sum(constructed.*G,2);


x = [ x_gap(1:i_gapStart-1) ; V ; x_gap(i_gapStart+gapLength:end) ]; 


end

