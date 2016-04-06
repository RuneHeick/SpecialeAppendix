close all;

Fs = 8000;                   % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 0.05;             % seconds
n = (-(StopTime-dt)/2:dt:(StopTime-dt)/2)';     % seconds

Fc = 300;                     % hertz
x = (sin(2*pi*linspace(100,Fc,size(n,1))'.*n))./n;
% x = (testdata(1431:1642))';
% x = x - mean(x);
n = (1:length(x))';

gapStart = 100;
gapSize = 20; 

 plot(n,x)
 hold on
 plot(n(gapStart:gapStart+gapSize),x(gapStart:gapStart+gapSize))

%% EMP Gap Splitting 


Left = [x(1:gapStart-1) , n(1:gapStart-1)]; 
Rigth = [ x(gapStart+gapSize:end) , n(gapStart+gapSize:end)];
signal = [Left ; Rigth];

infs = Modemd(x, gapStart, gapSize);
imfs2 = emd(x);

figure
index = 0; 
for i = 1:2:size(infs,1)*2
   index = index+1; 
   subplot(size(infs,1),2,i);
   plot(n, infs(index,:));
   hold on 
   plot(n(gapStart:gapStart+gapSize-1), infs(index,gapStart:gapStart+gapSize-1));
   
   if(size(imfs2,1)>= i)
       subplot(size(infs,1),2,i+1);
       plot(n, imfs2(index,:));
   end
end

figure
subplot(3,1,1);
plot(n,[ signal(1:length(Left),1) ; zeros(gapSize,1) ; signal(length(Left)+1:end,1) ])
subplot(3,1,2);
plot(n, x); 
subplot(3,1,3);
plot(n, sum(infs,1)); 