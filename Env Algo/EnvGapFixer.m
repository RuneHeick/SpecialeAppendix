function [ x_gap ] = EnvGapFixer( x_gap, i_gapStart, gapLength, traindelta )
%ENVGAPFIXER Summary of this function goes here
%   Detailed explanation goes here

n = (1:length(x_gap))';

i_trainLeft = (i_gapStart>traindelta)*(i_gapStart-traindelta) + (i_gapStart<=traindelta);
i_trainRigth = (length(x_gap)- (i_gapStart+gapLength) > traindelta)*(i_gapStart+gapLength+traindelta) + (length(x_gap)- (i_gapStart+gapLength) <= traindelta)*length(x_gap);

left = [x_gap(i_trainLeft:i_gapStart-1) , n(i_trainLeft:i_gapStart-1)]; 
rigth = [ x_gap(i_gapStart+gapLength:i_trainRigth) , n(i_gapStart+gapLength:i_trainRigth)];
signal = [left ; rigth];

%%
mins = []; 
maxes = [];

h = signal(:,1);
N = length(h);
      
[mins(2,:), mins(1,:)] = findpeaks(h*-1); 
 mins(2,:) = mins(2,:)*-1; 
[maxes(2,:), maxes(1,:)] = findpeaks(h); 

%Max
Lindex = find(maxes(1,:)<length(left)); 
Rindex = find(maxes(1,:)>length(left)+1);

if(isempty(Lindex) || isempty(Rindex))
  maxes = [ [1; h(1)] maxes [N; h(N)]];
else
maxesL = [[1 ; max(h(1),maxes(2,1))] maxes(:,Lindex) [length(left) ; max(h(length(left)),maxes(2,Lindex(end)))] ];
maxesR = [[length(left)+1 ; max(h(length(left)+1),maxes(2,Rindex(1)))]  maxes(:,Rindex) [N ; max(h(N),maxes(2,Rindex(end))) ] ];
maxes = [ maxesL maxesR];
end 

%Min 
Lindex = find(mins(1,:)<length(left)); 
Rindex = find(mins(1,:)>length(left)+1);
if(isempty(Lindex) || isempty(Rindex))
  mins = [ [1; h(1)] mins [N; h(N)]];
else
  minsL = [[1 ; min(h(1),mins(2,1))] mins(:,Lindex) [length(left) ; min(h(length(left)),mins(2,Lindex(end)))]];
  minsR = [[length(left)+1 ; min(h(length(left)+1),mins(2,Rindex(1)))]  mins(:,Rindex) [N ; min(h(N),mins(2,Rindex(end))) ] ];
  mins = [ minsL minsR];
end

%-------------------------------------------------------------------------
% spline interpolate to get max and min envelopes; form imf
maxenv = spline(signal(maxes(1,:),2),maxes(2,:),n(i_trainLeft:i_trainRigth));
minenv = spline(signal(mins(1,:),2), mins(2,:),n(i_trainLeft:i_trainRigth));

% figure(2)
% plot(n(i_trainLeft:i_trainRigth),maxenv)
% hold on 
% plot(n(i_trainLeft:i_trainRigth),minenv)
% plot(left(:,2), left(:,1))
% plot(rigth(:,2), rigth(:,1))
% hold off



x_gap(i_trainLeft:i_trainRigth) = imfRec2( x_gap(i_trainLeft:i_trainRigth)', i_gapStart-i_trainLeft, gapLength, minenv', maxenv' );

end

