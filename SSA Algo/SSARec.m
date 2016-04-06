


Fs = 8000;                   % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 0.05;             % seconds
n = (-(StopTime-dt)/2:dt:(StopTime-dt)/2)';     % seconds

Fc = 400;                     % hertz
x = (sin(2*pi*Fc*n)./n);
% x = (testdata(1431:1642))';
x = x - mean(x);
n = (1:length(x))';

gapStart = 100;
gapSize = 15; 

s = [ x(1:gapStart-1) ; zeros(gapSize,1) ; x(gapStart+gapSize:end) ];


%%
r2 = [];

Aerror = s(1:20); 

s = [ zeros(20,1); s(21:end) ]
best = 9999999999999999999999;
retSig = [];
sig = s; 
for i = 1:100

    while(1)
        sigold = sig; 
        [E,V,A,R] = ssa(s, 100); 

        sig = sum(R(:,1:i),2);
        
        figure(2)
        plot(sig)
        
        sig = sig-mean(sig);
        
        a = sum(abs(R(:,1:i)),2);
        b = sum(abs(R(:,i+1:end)),2);
        t = a+b; 

        ap = a./t; 
        bp = b./t; 
        
        sig = sig.*mean((1+bp));
          
        sig = [sig(1:20) ; s(21:gapStart-1) ; sig(gapStart:gapStart+gapSize-1) ; s(gapStart+gapSize:end) ];
       
        if(sum(abs(sig-sigold)) < 1e-5)
           break;  
        end
        
%         figure(1)
%         for t = 1:10
%             
%             subplot(10,1,t)
%             plot(R(:,t))
%         end
        
    end
    
    err = immse(sig(1:20),Aerror);
    
    if(err<best)
        retSig = sig + repmat(mean(s) - mean(sig), size(sig));
        best = err;
    else 
%         break; 
    end
    
end


figure(2)
plot(retSig);
hold on
plot(x, 'g')
hold off