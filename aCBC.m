function [C_ch1, C_ch2] = aCBC(red, green, Rmax, dr)

% INPUT (all valuesare in pixels or nm):
% - red: xy localisation coordinates of acquisition channel 1
% - green: xy localisation coordinates of acquisition channel 2
% - Rmax: maximal radius, in which colocalization will be probed. Value to be set
%         to the approximately the size of the expected structures/clusters
% - dr: binning radius for the distances distribution histograms. Value to
%       be set approximately to the resolution of the PALM/STORM image.


% OUTPUT:
% - C_ch1 = list of aCBC coefficients of localisations in acquisition channel 1
% - C_ch2 = list of aCBC coefficients of localisations in acquisition
% channel 2

r=linspace(.001,Rmax,round(Rmax/dr));

for ch = 1:2 % acquisition channel
    
    if ch == 1
        A = red;
        B = green;
    else
        A = green;
        B = red;
    end
    
    % Calculate distances (in px) between a localisation (Ap) and all
    % localisations from channel A & all localisations from channel B
    
    w = waitbar(0,['Calculating CBC values for channel ', num2str(ch),...
        ' of ', num2str(2),'...']);
    
    sizeA = size(A,2);
    E_ApB = zeros(1,sizeA);
    N_ApA_Rmax = zeros(1,sizeA);
    N_ApB_Rmax = zeros(1,sizeA);
    
    rho = zeros(1,sizeA);
    %     colocCoef = zeros(sizeA,1);
    
    for p = 1:sizeA
        
        waitbar(p/sizeA);
        
        % Measurement of distances between localisations:
        
        r_ApA = sqrt((A(1,p)-A(1,:)).^2 + (A(2,p)-A(2,:)).^2);
        r_ApA = [r_ApA(1:p-1) r_ApA(p+1:end)]; % removes the distance measurement of Ap with itself
        
        r_ApB = sqrt((A(1,p)-B(1,:)).^2 + (A(2,p)-B(2,:)).^2);
        
        
        % Distance b/n Ap and the nearest neighbour from species B :
        E_ApB(p) = min(r_ApB);
        
        % Number of localisations at Rmax:
        N_ApA_Rmax(p) = sum(r_ApA < Rmax);
        N_ApB_Rmax(p) = sum(r_ApB < Rmax);
        
        
        %         % Colocalisation coefficient :
        %         colocCoef(p) = sum(r_ApB < dr);
        
        % Number of localisations for each r:
        N_ApA = histc(r_ApA,r);
        N_ApB = histc(r_ApB,r);
        
        
        % Distribution of localisations for each r:
        D_ApA = (N_ApA(:)./(N_ApA_Rmax(p))).*(Rmax^2./r(:).^2);
        D_ApB = (N_ApB(:)./(N_ApB_Rmax(p))).*(Rmax^2./r(:).^2);
        
        
        % Rank of each distribution :
        [O_ApA, TIEADJ_ApA] = tiedrank(D_ApA);
        [O_ApB, TIEADJ_ApB] = tiedrank(D_ApB);
        
        % Spearman's correlation coefficient :
        [rho_p, pval] = corr([O_ApA O_ApB], 'type', 'Spearman');
        
        if isnan(rho_p(1,2))
            rho(p) = 0;
        else
            rho(p) = rho_p(1,2);
        end        
    end
    close(w);
    
    % Colocalisation values for each particle of species A :
    C_A = rho.*exp(- E_ApB./Rmax) ;
    
    if ch == 1
        C_ch1 = C_A; % Colocalisation values for particles from the red channel
        %         ColocCoef_ch1 = colocCoef./max(colocCoef);
    else
        C_ch2 = C_A; % Colocalisation values for particles from the green channel
        %         ColocCoef_ch2 = colocCoef./max(colocCoef);
    end    
end
end