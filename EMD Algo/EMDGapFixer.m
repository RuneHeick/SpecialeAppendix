function [ x_gap ] = EMDGapFixer( x_gap, i_gapStart, gapLength, traindelta )
%EMDGAPFIXER Summary of this function goes here
%   Detailed explanation goes here

n = (1:length(x_gap))';

i_trainLeft = (i_gapStart>traindelta)*(i_gapStart-traindelta) + (i_gapStart<=traindelta);
i_trainRigth = (length(x_gap)- (i_gapStart+gapLength) > traindelta)*(i_gapStart+gapLength+traindelta) + (length(x_gap)- (i_gapStart+gapLength) <= traindelta)*length(x_gap);

i_gapZoomStart = i_gapStart - i_trainLeft +1;
x_gapZoom = x_gap(i_trainLeft:i_trainRigth);
%%

x_gapZoomMean = mean(x_gapZoom); 
x_gapZoom = x_gapZoom - x_gapZoomMean;

%%

%returns fixed imfs 
imfs = Modemd(x_gapZoom, i_gapZoomStart, gapLength);
recoverdSignal = sum(imfs,1); 

x_gap(i_gapStart:i_gapStart+gapLength-1) = recoverdSignal(i_gapZoomStart:i_gapZoomStart+gapLength-1)+x_gapZoomMean; 


end

