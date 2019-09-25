function [iv, ki]=OM_Interpolate_IV(mnes, vol, gridparams)
% version 2012-08-07
% Please write to vilkov@vilkov.net with questions/comments
% 2014 Grigory Vilkov


% set the grid in terms of moneyness that we want to use to compute the
    % MFIV/ MFIS and other integrals as needed
    m = gridparams.m; % use points -500:500, i.e. 1001 points for integral approximation
    k = gridparams.k; % use moneyness 1/(1+2) to 1+2 - should be enough as claimed by JT
    u = (1+k)^(1/m);
    %create strikes using s=1, i.e. get k/s = moneyness
    mi = -m:m;
    ki = u.^mi';
    iv = zeros(size(ki))*NaN; % preallocate the implied volatilities vector


    currspline = pchip(mnes,vol);

    % now get the inter-/ extrapolated BS IV for the set of strikes

    k_s_max = max(mnes); % OTM calls
    k_s_min = min(mnes); % OTM puts

    iv_max=vol(1); %for more OTM puts i.e we have iv_max for min k/s, i.e. for OTM put option
    iv_min=vol(end);%for more OTM calls  i.e. we have iv_min for OTM call option

    %calculate the interpolated/ extrapolated IV for these k/s
    ks_larger_ind   = (ki>k_s_max); % more OTM call
    ks_smaller_ind  = (ki<k_s_min); % more OTM put
    ks_between_ind  = (ki>=k_s_min & ki<=k_s_max);
    
    if ~isempty(ks_larger_ind)
        iv(ks_larger_ind,1)=iv_min;  % more OTM put get iv_min
    end
    if ~isempty(ks_smaller_ind)
        iv(ks_smaller_ind,1)=iv_max; % more OTM call get iv_max
    end
    if ~isempty(ks_between_ind)
        iv(ks_between_ind,1)=ppval(currspline,ki(ks_between_ind));
    end