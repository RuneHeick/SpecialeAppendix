

close all;

Fs = 8000;                   % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 0.05;             % seconds
n = (-(StopTime-dt)/2:dt:(StopTime-dt)/2)';     % seconds

Fc = 200;                     % hertz
x = (sin(2*pi*linspace(100,Fc,size(n,1))'.*n));
% x = (testdata(1431:1642))';
% x = x - mean(x);
% n = (1:length(x))';

gapStart = 105;
gapSize = 25; 


%% 

left = [x(1:gapStart-1) , n(1:gapStart-1)]; 
rigth = [ x(gapStart+gapSize:end) , n(gapStart+gapSize:end)];
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
maxenv = spline(signal(maxes(1,:),2),maxes(2,:),n);
minenv = spline(signal(mins(1,:),2), mins(2,:),n);

xnew = imfRec2( x', gapStart, gapSize, minenv', maxenv' )

figure(1)
plot(n,xnew);
hold on 
plot(n,x);
plot(n,maxenv);
plot(n,minenv); 
hold off;
