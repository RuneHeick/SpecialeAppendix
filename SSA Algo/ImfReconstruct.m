function [ inf ] = ImfReconstruct( inf, gapStart, gapLength, envlower, envUpper )
%INFRECONSTRUCT Summary of this function goes here
%   Detailed explanation goes here

%% Find I 

    
    [vmax, maxs ] = findpeaks(inf);
    [vmin, mins ] = findpeaks(inf*-1);
    
    if(max(vmax)<0.001 && min(vmin)>-0.001)
        return; 
    end
    
    
    t_left = mins(max(find(mins<(gapStart-1)))); 
    T_left = maxs(max(find(maxs<(gapStart-1))));
    
    t_rigth = mins(min(find(mins>gapStart+gapLength))); 
    T_rigth = maxs(min(find(maxs>gapStart+gapLength)));
    
    if(length(maxs)<4 || isempty(t_left) || isempty(T_left) || isempty(t_rigth) || isempty(T_rigth) )
        %Do other reconstruction 
        ipoint = [ gapStart-2 gapStart-1 gapStart+gapLength gapStart+gapLength+1]; 
        vpoint = inf(ipoint);
        
        inf(gapStart:gapStart+gapLength) = spline(ipoint,vpoint,gapStart:gapStart+gapLength);
        return; 
    end

    
    a_left = ((inf(T_left)-inf(t_left))/(T_left-t_left));
    k_left = -floor((inf(t_left) - a_left*t_left)/a_left);

    a_rigth = ((inf(T_rigth)-inf(t_rigth))/(T_rigth-t_rigth));
    k_rigth = -floor((inf(t_rigth) - a_rigth*t_rigth)/a_rigth);

    I = [];

    I1 = max([t_left, T_left, k_left]); 
    I = [I I1]; 

    if I1 == k_left && k_left>T_left && k_left>t_left
       I2 = max(t_left, T_left);
       I = [I2 I]; 
    end

    I3 = min([t_rigth,T_rigth,k_rigth]);
    I = [I I3]; 

    if I3 == k_rigth && k_rigth<t_rigth && k_rigth<T_rigth
       I4 = min(t_rigth,T_rigth);
       I = [I I4]; 
    end

    %% Find T

    P_left = 1;
    P_rigth = 1;
    
    q = gapStart;
    Q = gapLength;
    
    d_left = abs(t_left-T_left); 
    d_rigth = abs(t_rigth-T_rigth); 
    
    d = min(d_left,d_rigth);

    D_left = q - max(T_left,t_left); 
    D_rigth = min(T_rigth,t_rigth)-(q+Q-1); 
    
    J_left = q - floor(P_left*D_left); 
    J_rigth = q+Q-1+floor(P_rigth*D_rigth); 
    
    Tdelta = floor((J_rigth - J_left)/d);
    TMark = Tdelta-1;
  
    %% Create I-Lib 
    I_lib = {...
    [t_left, t_rigth ] [1 0] [1 1] [2 1]
    [t_left, k_rigth, t_rigth] [1 0] [1 1] [2 1]
    [t_left, k_left,k_rigth, t_rigth ] [1 0] [1 1] [2 1]
    [t_left, k_left, t_rigth ] [1 0] [1 1] [2 1]
    %Next bin 
    [t_left, T_rigth] [1 0] [1 0] [2 0]
    [t_left, k_rigth,T_rigth] [1 0] [1 0] [2 0]
    [T_left, k_left, T_rigth] [1 0] [1 0] [2 0]
    [t_left, k_left, T_rigth] [1 0] [1 0] [2 0]
    [T_left, k_rigth,t_rigth] [1 0] [1 0] [2 0]
    [T_left, t_rigth] [1 0] [1 0] [2 0]
    %Next bin 
    [T_left, T_rigth] [1 1] [1 0] [2 1]
    [T_left, k_left, k_rigth, T_rigth] [1 1] [1 0] [2 1]
    [T_left, k_rigth, T_rigth] [1 1] [1 0] [2 1]
    [T_left, k_left, T_rigth] [1 1] [1 0] [2 1]
    %Next bin 
    [T_left, k_left, k_rigth, t_rigth] [1 1] [1 1] [2 2]
    [t_left, k_left, k_rigth, T_rigth] [1 1] [1 1] [2 2]
    };
    
%% Find I in I_lib
    
    LibIndex = find(cellfun(@(x)(isequal(x,I)), I_lib(:,1))); 
    ABT = I_lib(LibIndex,2:end); 

    h = ceil((TMark-ABT{1,3}(2))/(ABT{1,3}(1))); 
    
    T_Hat = ABT{1,3}(2)+ABT{1,3}(1)*h; 
    k = 1:ceil(T_Hat); 
    indexs = J_left+k*mean(d_left,d_rigth);;
    
    if(~isempty(indexs))
        i_val = (indexs>gapStart-1) .* (indexs<gapStart+gapLength);
        indexs = indexs(find(i_val));

        SetMin = T_left>t_left;
        if(~isempty(indexs))
            for i = 1:size(indexs,2)
                if SetMin
                    inf(indexs(i)) = envlower(indexs(i)); 
                    SetMin = 0; 
                else
                    inf(indexs(i)) = envUpper(indexs(i)); 
                    SetMin = 1; 
                end
            end
        end
    end
    ipoint = unique([t_left T_left gapStart-1 indexs gapStart+gapLength t_rigth, T_rigth],'sorted'); 
   
    vpoint = inf(ipoint);
    
    figure(7)
    plot(inf)
    hold on 
    plot(envlower)
    plot(envUpper)
    scatter(ipoint,vpoint);
    
    
    
    inf(gapStart:gapStart+gapLength) = spline(ipoint,vpoint,gapStart:gapStart+gapLength);
    
    plot(inf)
    hold off
    
end

