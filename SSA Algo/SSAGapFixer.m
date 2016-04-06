function [ x_gap ] = SSAGapFixer( x_gap, i_gapStart, gapLength, traindelta )
%SSAGAPFIXER Summary of this function goes here
%   Detailed explanation goes here

i_trainLeft = (i_gapStart>traindelta)*(i_gapStart-traindelta) + (i_gapStart<=traindelta);
i_trainRigth = (length(x_gap)- (i_gapStart+gapLength) > traindelta)*(i_gapStart+gapLength+traindelta) + (length(x_gap)- (i_gapStart+gapLength) <= traindelta)*length(x_gap);


knownSignal = [ x_gap(i_trainLeft:i_gapStart-1) ; x_gap(i_gapStart+gapLength:i_trainRigth) ];

signalMean = mean(knownSignal); 

s = [ x_gap(i_trainLeft:i_gapStart-2)-signalMean ; zeros(gapLength,1) ; x_gap(i_gapStart+gapLength+1:i_trainRigth)-signalMean ];
i_gapstartS = i_gapStart-i_trainLeft;

%%

AerrorLen = min(20, floor(traindelta/4));

if(AerrorLen == 0)
   fix = interp1q([1 2+gapLength]',[x_gap(i_gapStart-1) x_gap(i_gapStart+gapLength)]',[2:1+gapLength]');
   x_gap = [ x_gap(1:i_gapStart-1); fix ; x_gap(i_gapStart+gapLength:end) ];
   return;
end

Aerror = s(1:AerrorLen); 
s = [ zeros(AerrorLen,1); s(AerrorLen+1:end) ];

best = 9999999999999999999999;
retSig = [];
sig = s; 
lag = min(100,traindelta); 
for i = 1:lag

    while(1)
        sigold = sig; 
        [E,V,A,R] = ssa(s, lag); 

        sig = sum(R(:,1:i),2);
        
        sig = sig-mean(sig);
        
        a = sum(abs(R(:,1:i)),2);
        b = sum(abs(R(:,i+1:end)),2);
        t = a+b; 

        if t ~= 0
            ap = a./t; 
            bp = b./t; 

            sig = sig.*mean((1+bp));
        end
        
        sig = [sig(1:AerrorLen) ; s(AerrorLen+1:i_gapstartS-2) ; sig(i_gapstartS-1:i_gapstartS+gapLength) ; s(i_gapstartS+gapLength+1:end) ];
       
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
    
    err = immse(sig(1:AerrorLen),Aerror);
    
    if(err<best)
        retSig = sig;
        best = err;
    else 
         break; 
    end
    
end

offset = mean([retSig(i_gapstartS-1)-x_gap(i_gapStart-1), retSig(i_gapstartS+gapLength)-x_gap(i_gapStart+gapLength)]);
 
x_gap(i_gapStart:i_gapStart+gapLength-1) = retSig(i_gapstartS:i_gapstartS+gapLength-1)-offset;

end

