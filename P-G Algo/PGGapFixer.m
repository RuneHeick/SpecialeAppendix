function [ x_gap ] = PGGapFixer( x_gap, i_gapStart, gapLength, traindelta )
%PGGAPFIXER Summary of this function goes here
%   Detailed explanation goes here

i_trainLeft = (i_gapStart>traindelta)*(i_gapStart-traindelta) + (i_gapStart<=traindelta);
i_trainRigth = (length(x_gap)- (i_gapStart+gapLength) > traindelta)*(i_gapStart+gapLength+traindelta) + (length(x_gap)- (i_gapStart+gapLength) <= traindelta)*length(x_gap);

i_gapZoomStart = i_gapStart - i_trainLeft +1;
x_gapZoom = x_gap(i_trainLeft:i_trainRigth);
x_gapZoomMean = mean(x_gapZoom); 
x_gapZoom = x_gapZoom -x_gapZoomMean;

D = diag([ones(i_gapZoomStart-1,1) ; zeros(gapLength,1) ; ones(size(x_gapZoom,1)-(i_gapZoomStart-1)-gapLength,1)]);

y = D*x_gapZoom; 
bestfit = 0; 
isStarted = 0; 
bestfitpoint = 0; 

% figure
% plot(abs(fft(y)))


for M = 1:length(y)/2
    gamma = [ones(M,1) ; zeros(size(x_gapZoom,1)-2*M,1) ; ones(M,1)];
    GAMMA = diag(gamma);

    F = dftmtx(size(x_gapZoom,1)); 
    B = inv(F)*GAMMA*F; 
    I = eye(size(x_gapZoom,1));


    x_hat = y;
    x_old = y; 

    delta = 0.000001; 

        for i = 1:100

            x_hat = y + (I-D)*B*x_hat;

            if(sum(abs(x_old-x_hat))/size(x_hat,1) < delta)
                besthat = x_hat;
                break;  
            end

            x_old = x_hat; 

        end

    sigFFT = fft(real(x_hat));

    E_t = sum(abs(sigFFT(1:end/2)).^2);
    E_h = sum(abs(sigFFT(ceil(end/2-end/5):end/2)).^2);
 
    
    g(M) = log10(E_h/E_t);
    
%     if(M>2 && abs(g(M))-abs(gold) > bestfit )
%         bestfit = abs(g(M))-abs(gold);
%         besthat = x_hat;
%         bestfitpoint = g(M);
%         isStarted = 1; %true
%         bestM = M; 
%     elseif(isStarted && g(M) < bestfitpoint )
%         besthat = x_hat;
%         bestfitpoint = g(M);
%         bestM = M;
%     elseif(isStarted && g(M) > bestfitpoint )
%         isStarted = 0; 
    if(M>2 && g(M-1) < bestfit )
        bestfit = g(M-1);
        besthat = x_hat;
        bestM = M-1; 
    end
    gold = g(M);
end
% 
% figure(2)
% plot(g)

x_gap(i_trainLeft:i_trainRigth) = real(besthat)+x_gapZoomMean;


end

