function [ imf ] = imfRec2( imf, gapStart, gapLength, envlower, envUpper )
%IMFREC2 Summary of this function goes here
%   Detailed explanation goes here


    [vmax, maxs ] = findpeaks(imf);
    [vmin, mins ] = findpeaks(imf*-1);
     
    t_left = mins(find(mins<(gapStart-1)));
    if(imf(gapStart-1)<envlower(gapStart-1) )
        t_left = [ t_left gapStart-1];
    end
    t_left(2,:) = zeros(1,size(t_left,2));
    
    T_left = maxs(find(maxs<(gapStart-1)));
    if(imf(gapStart-1)>envUpper(gapStart-1))
        T_left = [ T_left gapStart-1];
    end
    T_left(2,:) = ones(1,size(T_left,2));
    
    t_rigth = mins(find(mins>gapStart+gapLength)); 
    if(imf(gapStart+gapLength)<envlower(gapStart+gapLength) )
        t_rigth = [ gapStart+gapLength t_rigth];
    end
    t_rigth(2,:) = zeros(1,size(t_rigth,2));
    
    T_rigth = maxs(find(maxs>gapStart+gapLength));
    if(imf(gapStart+gapLength)>envUpper(gapStart+gapLength) )
        T_rigth = [gapStart+gapLength T_rigth];
    end
    T_rigth(2,:) = ones(1,size(T_rigth,2));
    
    if( size(t_left,2)+size(T_left,2) < 2 || size(t_rigth,2)+size(T_rigth,2) < 2 )
        %Do other reconstruction 
        ipoint = [ gapStart-2 gapStart-1 gapStart+gapLength gapStart+gapLength+1]; 
        vpoint = imf(ipoint);

        imf(gapStart:gapStart+gapLength) = spline(ipoint,vpoint,gapStart:gapStart+gapLength);
        return; 
    end

    leftPeaks = sortrows ([t_left T_left]')'; 
    rigthPeaks = sortrows([t_rigth T_rigth]')'; 

    peakmap = [ leftPeaks(1,2:end) - leftPeaks(1,1:end-1) rigthPeaks(1,2:end) - rigthPeaks(1,1:end-1); 
                leftPeaks(1,1:end-1) rigthPeaks(1,1:end-1);
                leftPeaks(2,1:end-1) rigthPeaks(2,1:end-1)
                ];
        
    p = polyfit(peakmap(2,:),peakmap(1,:),max(ceil(size(peakmap,2)-2),1));
    y1 = polyval(p,peakmap(2,:));
    
%     figure(4)
%     plot(peakmap(2,:),y1);
%     hold on 
%     plot(peakmap(2,:),peakmap(1,:));
%     hold off
    
    error = abs(peakmap(1,:)-y1);
    div = std(error);
       
    peak = leftPeaks(1,end);
    rconp = [];
    while(peak < gapStart+gapLength)
        if(peak>gapStart)
            if(isempty(rconp))
               rconp = [peak ; ~leftPeaks(2,end)];
            else
                rconp = [rconp [peak ; ~rconp(2,end)]];
            end            
        end
        
        prediction = polyval(p,peak);
        
        prediction = (prediction>=1)*prediction + (prediction<1);
        peak = floor((peak + prediction) + div.*randn(1,1));
    end
    
    if(isempty(rconp))
        ipoint = unique([ t_left(1,:) gapStart-2 gapStart-1 T_left(1,:) gapStart+gapLength gapStart+gapLength+1 T_rigth(1,:) ], 'sorted');
    else
        
        if( (leftPeaks(2,end) == rigthPeaks(2,1) && mod(size(rconp,2),2) == 0) || (leftPeaks(2,end) ~= rigthPeaks(2,1) && mod(size(rconp,2),2) == 1) )
        rconp = rconp(:,1:end-1); 
        end
    
        upperIndex = rconp(1,find(rconp(2,:)==1)); 
        lowerIndex = rconp(1,find(rconp(2,:)==0));

        imf(lowerIndex) = envlower(lowerIndex); 
        imf(upperIndex) = envUpper(upperIndex); 

        ipoint = unique([t_left(1,:) gapStart-2 gapStart-1 lowerIndex t_rigth(1,:) T_left(1,:) gapStart+gapLength gapStart+gapLength+1 upperIndex T_rigth(1,:) ], 'sorted'); 
        
    end
    
    vpoint = imf(ipoint);

%     figure(2)
%     plot(imf)
%     hold on
%     plot(envlower)
%     plot(envUpper)
    
    imf(gapStart:gapStart+gapLength) = spline(ipoint,vpoint,gapStart:gapStart+gapLength);
    
%     plot(imf)
%     hold off;
    
end

