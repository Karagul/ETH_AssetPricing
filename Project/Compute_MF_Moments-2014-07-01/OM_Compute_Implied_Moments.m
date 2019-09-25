function [moments, nopt] = OM_Compute_Implied_Moments(optdata, rate, gridparams)
% ####################################################
% COMPUTE IMPLIED MOMENTS FUNCTION: 
% USING OWN ROUTINE WITH INTERP: MFIV(variance)/MFIS/MFIK/CXI(variance) (corridor MFIV with interpolation)
% USING CBOE ROUTINE:            VIX(as variance)/ SKEW
% version July 01, 2014
% Please write to vilkov@vilkov.net with questions/comments
% 2014 Grigory Vilkov


% ## OPTIONS: INPUT FILE
% 1 b.secid,
% 2 b.date,
% 3 b.exdate-b.date as days,
% 4 b.last_date,
% 5 b.strike_price/1000/abs(b1.close) as mnes,
% 6 b.best_bid/abs(b1.close),
% 7 b_best_offer/abs(b1.close),
% 8 b.volume,
% 9 b.open_interest,
% 10 b.impl_volatility,
% 11 b.delta
% 
% OIMeasureNames  = {
%     'mfiv'          %1 mfiv with interp              (Bakshi/Kapadia/Madan)
%     'smfiv'         %2 simple var swap^2 with interp (Martin 2012)
%     'mfis'          %3 mfis with interp              (Bakshi/Kapadia/Madan)
%     'mfik'          %4 mfik with interp              (Bakshi/Kapadia/Madan)
%     'cxi'           %5 corridor mfiv with interp     (Andersen/Bondarenko)
%     'tlmi'          %6 tail loss measure with interp (Hamidieh)
%     'vix'           %7 vix^2: not interpolated       (CBOE white paper)
%     'svix'          %8 simple var swap^2 not interp  (Martin 2012)
%     'skew'          %9 skew: just S, not 100-10S yet (CBOE white paper)
%     'cx'            %10 corridor mfiv w/o interp     (Andersen/Bondarenko)
%     'tlm'           %11 tail loss measure w/o interp (Hamidieh)
%     'iv'   };       %12 ATM IV  
% 

options = optimset('TolFun', 1e-16, 'TolX', 1e-16,'MaxIter',5000,'MaxFunEvals',100000,'Display','Off');
% default params
if nargin<3
    gridparams.m            = 500;
    gridparams.k            = 2;
    gridparams.meascompute  = [1:5 7 8 9 10 12];
    gridparams.filter       = 2;
    gridparams.mnes_limits  = [0.7 1.3];        % filter out the options beyond these mnes limits, by default [0.7 1.3]
    gridparams.cx_limits    = [0.03 0.97];
    gridparams.tlm_limits   = 0.4; % delta limit for the tail measure

OIMeasureNames  = {
    'mfiv'          %1 mfiv with interp              (Bakshi/Kapadia/Madan)
    'smfiv'         %2 simple var swap^2 with interp (Martin 2012)
    'mfis'          %3 mfis with interp              (Bakshi/Kapadia/Madan)
    'mfik'          %4 mfik with interp              (Bakshi/Kapadia/Madan)
    'cxi'           %5 corridor mfiv with interp     (Andersen/Bondarenko)
    'tlmi'          %6 tail loss measure with interp (Hamidieh)
    'vix'           %7 vix^2: not interpolated       (CBOE white paper)
    'svix'          %8 simple var swap^2 not interp  (Martin 2012)
    'skew'          %9 skew: just S, not 100-10S yet (CBOE white paper)
    'cx'            %10 corridor mfiv w/o interp     (Andersen/Bondarenko)
    'tlm'           %11 tail loss measure w/o interp (Hamidieh)
    'iv'   };       %12 ATM IV                         ()
     
    gridparams.measnames    = OIMeasureNames(:,1);
    gridparams.surface      = false;
end

% predefine the output to be NaN
if nargout==2
    nopt=NaN;
end

for z=1:length(gridparams.measnames)
    moments.(gridparams.measnames{z}) = NaN;
end

% check what routine to use CBOE or OURs
ComputeOURStuff  = ~isempty(intersect(gridparams.meascompute,[1:6]));
ComputeCBOEStuff = ~isempty(intersect(gridparams.meascompute,[7:12]));


% ##############################
mat      = optdata(1,3)/365; % days to maturity
er       = exp(mat*rate);

if ComputeOURStuff    
    % APPLY THE FILTER FOR INTERPOLATION METHODS
    % ### define and apply filters
    
    % FILTER 1: REMOVE
    %         zero bids
    %         mnes<lim1 | mnes > lim2
    %         OTM only - by delta
    
    % FILTER 2: REMOVE
    %         zero bids
    %         zero open interest
    %         mnes<lim1 | mnes > lim2
    %         OTM only - by delta
    
    % FILTER 3: REMOVE
    %         zero bids
    %         zero open interest
    %         mnes<lim1 | mnes > lim2
    %         OTM only - by delta
    %         prc<0.5
    
    lim1 = gridparams.mnes_limits(1);        % filter out the options beyond these mnes limits, by default [0.7 1.3]
    lim2 = gridparams.mnes_limits(2);
    switch gridparams.filter
        case 1
            bads= optdata(:,6)==0 |                   optdata(:,5) <lim1 | optdata(:,5)>lim2 | optdata(:,11)>=0.5 | optdata(:,11)<-0.5 ;
        case 2
            bads= optdata(:,6)==0 | optdata(:,9)==0 | optdata(:,5) <lim1 | optdata(:,5)>lim2 | optdata(:,11)>=0.5 | optdata(:,11)<-0.5;
        case 3
            bads= optdata(:,6)==0 | optdata(:,9)==0 | optdata(:,5) <lim1 | optdata(:,5)>lim2 | (optdata(:,6)+optdata(:,7))/2<0.5 | optdata(:,11)>=0.5 | optdata(:,11)<-0.5;
    end
    
    % remove the bad observations
    optdata4interp = optdata(~bads,:);
    
    % moneyness is sorted now from low to high
    optdata4interp  = sortrows(optdata4interp,5);
    
    
    % take the first from any duplicate points
    dupl = optdata4interp(1:end-1,5) == optdata4interp(2:end,5);
    dupl = [false;dupl]; 
    optdata4interp(dupl,:) = [];
        
    
    mnes            = optdata4interp(:,5);
    vol             = optdata4interp(:,10);
    nopt            = length(mnes);
    % now check how many options we have, and if fewer than 4, do not compute anyting
    % of there are duplicate moneyness points, do not compute anything
    if length(optdata4interp(:,1))>=4 && length(unique(optdata4interp(:,5)))==length(optdata4interp(:,5))% then proceed
        
        % ######################################
        % compute the moments WITH INTERPOLATION
        % ######################################
        [iv, ki] = OM_Interpolate_IV(mnes, vol, gridparams);
        
        % calculate a bunch of BS prices
        % ###         OTM calls
        otmcalls = ki>=1;
        [currcalls, itmputs] = blsprice(1,ki(otmcalls),rate,mat,iv(otmcalls));
        % ###          OTM puts
        kiputs                = ki(~otmcalls);
        [itmcalls, currputs]  = blsprice(1,kiputs,rate,mat,iv(~otmcalls));
        [~, currputd]         = blsdelta(1,kiputs,rate,mat,iv(~otmcalls));
        
        % use the gridparams to define some variables
        m   = gridparams.m; % use points -500:500, i.e. 1001 points for integral approximation
        k   = gridparams.k; % use moneyness 1/(1+2) to 1+2 - should be enough as claimed by JT
        u   = (1+k)^(1/m);
        mi  = -m:m;
        
        % ### compute MFIV
        a   = 2*(u-1);
        ic  = mi(otmcalls)';
        ip  =(mi(~otmcalls))';
        b1  = sum((1-(log(1+k)/m).*ic).*currcalls./u.^ic);
        b2  = sum((1-(log(1+k)/m).*ip).*currputs./u.^ip);
        V   = a*(b1+b2);
        moments.mfiv=V/mat;
        
        % ### compute SMFIV
        Ki       = u.^[ic; ip]; % get all moneyness points        
        Q        = [currcalls; currputs];        
        dKi      = NaN(size(Q,1),1);
        dKi(2:end-1,1)    = (Ki(3:end) - Ki(1:end-2))/2;
        dKi(1)   = Ki(2) - Ki(1);
        dKi(end) = Ki(end) - Ki(end-1);
        
        K0sq    = 1; 
        
        % compute SVIX
        moments.smfiv = (2 * er * sum(dKi.*Q./K0sq))/mat;
        
        
        % ### compute MFIS/MFIK
        a   = 3*(u-1)*log(1+k)/m;
        ic  = mi(otmcalls)';
        ip  = (mi(~otmcalls))';
        b1  = sum(ic.*(2-(log(1+k)/m).*ic).*currcalls./u.^ic);
        b2  = sum(ip.*(2-(log(1+k)/m).*ip).*currputs./u.^ip);
        W   = a*(b1+b2);
        
        a   = 4*(u-1)*(log(1+k)/m)^2;
        ic  = mi(otmcalls)';
        ip  = (mi(~otmcalls))';
        b1  = sum(ic.^2.*(3-(log(1+k)/m).*ic).*currcalls./u.^ic);
        b2  = sum(ip.^2.*(3-(log(1+k)/m).*ip).*currputs./u.^ip);
        X   = a*(b1+b2);
        
        mu  = er-1-er/2*V-er/6*W-er/24*X;
        c   = (er*V-mu^2);
        
        moments.mfis = (er*W-3*mu*er*V+2*mu^3)/(c^(3/2)); % hope, it's correct :D:D:
        moments.mfik = (er*X-4*mu*er*W+6*er*mu^2*V-3*mu^4)/c^2;
        
            
        if sum(strcmpi('tlmi',gridparams.measnames(gridparams.meascompute)))>0
            % ### compute TLMI
            if gridparams.tlm_limits>0.5
                goods = kiputs <= 1-gridparams.tlm_limits*sqrt(moments.mfiv/251);
            else                
                goods = abs(currputd)<=gridparams.tlm_limits;
            end
            if sum(goods)>1
                Pt          = currputs(goods); % prices
                K           = kiputs(goods);   % strikes
                [P0,index]  = max(Pt);
                K0          = K(index);
                
                optResult   = fminsearch(@(x) objFunc_TailMeasure(x,P0, K0, Pt,K),[0.5,0.01], options);
                moments.tlmi= optResult(2)./(1-optResult(1)); % optResult(1);%
            end
        end
        
        % ### compute CX
        
        % identify the break point
        % R = P(K)/(P(K)+C(K)) => K for R = gridparams.cx_limits(1) and gridparams.cx_limits(2)
        
        R =  currputs./(currputs+itmcalls); goodputs  = R>=gridparams.cx_limits(1); currputs  = currputs(goodputs);
        R =  itmputs./(currcalls+itmputs);  goodcalls = R<=gridparams.cx_limits(2); currcalls = currcalls(goodcalls);
        
        % compute the MFIV -> CX using only the qualified calls/puts
        a   = 2*(u-1);
        ic  = mi(otmcalls)';    ic = ic(goodcalls);
        ip  = (mi(~otmcalls))'; ip = ip(goodputs);
        b1  = sum((1-(log(1+k)/m).*ic).*currcalls./u.^ic);
        b2  = sum((1-(log(1+k)/m).*ip).*currputs./u.^ip);
        V   = a*(b1+b2);
        moments.cxi=V/mat;
        
        
    end
end  % ComputeOURStuff

% if isnan(moments.mfiv); pause; end; % for debugging -- 

% #########################################
% compute the moments WITHOUT INTERPOLATION
% #########################################
if ComputeCBOEStuff
    % #### define the ATM level by the smallest call/put price diff
    % The at- the-money strike is defined as the listed strike immediately below
    % the S&P 500 forward price. To find the forward price, find the strike for which
    % the difference between the midquotes of the call and put is at a minimum. Then
    % calculate the forward price as F=erT *(C-P)+K,where T is the time to expiration, C and P are
    % the call and put midquotes, and K is the strike at which minimum occurs.
   
    
    % ## allocate calls and puts by mnes
    mnesall      = unique(optdata(:,5));
    opttemp      = NaN(length(mnesall),7);
    opttemp(:,1) = mnesall;
    flds         = {'pts','cls'};
    temp.pts     = sortrows(optdata(optdata(:,11)<0,:),5); % puts
    temp.cls     = sortrows(optdata(optdata(:,11)>0,:),5); % calls
    for fld = 1:2 % go thru puts and calls
        tmp = temp.(flds{fld}); % puts OR calls
        for z=1:size(tmp,1)
            curr = tmp(z,5) == mnesall; % select option with a given mnes level
            if sum(curr)==1 % and if there is only one option with this mnes, save mid-prices and bids and IV
                opttemp(curr,fld+1) = (tmp(z,6)+tmp(z,7))/2; % mid
                opttemp(curr,fld+3) = tmp(z,6);  % bid
                opttemp(curr,6)     = tmp(z,10); % iv
                opttemp(curr,7)     = tmp(z,11); % delta
            end
        end
    end
    
    % ## find the point with the min diff between calls and puts
    if gridparams.surface
        fut = 1/er;
    else
        diff        = (opttemp(:,3)-opttemp(:,2));
        [~, atmind] = min(abs(diff)); atmind = atmind(1);
        atmmnes     = opttemp(atmind,1);
        diff        = diff(atmind);
        % check how far the ATM strike is from 1
        % if more than 3%, use the mnes = 1/er to get the ATM option
        if abs(atmmnes-1)>0.03
            fut = 1/er;
        else
            fut = diff*er + atmmnes;
        end
    end
    
    % find ATM strike defined as the listed strike immediately below the S&P 500 forward price.
    atm = find(opttemp(:,1)<=fut,1,'last'); atm = opttemp(atm,1); if isempty(atm); atm=1/er; end;
    % use the atm price to find OTM options
    allputs  = opttemp(:,1)<=atm; allputs  = opttemp(allputs, [1 4 2 6 7]);  % mnes bid mid iv delta
    allcalls = opttemp(:,1)>=atm; allcalls = opttemp(allcalls,[1 5 3 6 7]);  % mnes bid mid iv delta
    allputs(isnan(allputs))   = 0; 
    allcalls(isnan(allcalls)) = 0;
    % ## find two consecutive zero bids and eliminate everythin beyond
    % sort put by mnes in ascending order, i.e.,
    % for puts:
    % 0.7 0.8 0.9 0.95 1.0... mnes: more OTM -> ATM
    allputs   = sortrows(allputs,1);   
    zerobids  = allputs(:,2)==0;
    bids2zero = find(zerobids(1:end-1)+zerobids(2:end) == 2,1,'last');     
    zerobids(1:bids2zero) = 1;
    allputs(zerobids,:)   = []; % eliminate all bad puts
    
    % sort calls by mnes in descending order, i.e.,
    % for calls:
    % 1.2 1.1 1.05 1.02 1.0... mnes: more OTM -> ATM
    allcalls  = sortrows(allcalls,-1);
    zerobids  = allcalls(:,2)==0;
    bids2zero = find(zerobids(1:end-1)+zerobids(2:end) == 2,1,'last'); 
    zerobids(1:bids2zero) = 1;
    allcalls(zerobids,:)  = []; % eliminate all bad calls
    allcalls   = sortrows(allcalls,1); % sort again in the ascending order
    
    % take options closest to ATM and compute average IV
    iv = NaN(2,1);
    if size(allputs,1)>0
        iv(1) = allputs(end,4);
    end
    if size(allcalls,1)>0
        iv(2) = allcalls(end,4);
    end
    % compute iv
    moments.iv = nanmean(iv.^2); if moments.iv==0; moments.iv = NaN; end; 
    
    % ### coompute if at least 2 call or 2 puts
    if size(allcalls,1)>1 && size(allputs,1)>1
        % ## prepare all the variables for computation
        %?K for the lowest strike is simply the difference between the lowest strike
        % and the next higher strike. Likewise, ?K for the highest strike is the difference
        % between the highest strike and the next lower strike.
        Q       = [allputs(:,[1 3]); allcalls(:,[1 3])]; % all options: sorted as OTM puts -> OTM calls
        Ki       = Q(:,1);
        Q        = Q(:,2);
        dKi      = NaN(size(Q,1),1);
        dKi(2:end-1,1)    = (Ki(3:end) - Ki(1:end-2))/2;
        dKi(1)   = Ki(2) - Ki(1);
        dKi(end) = Ki(end) - Ki(end-1);
        
        F0 = fut;
        K0 = atm;
        
        fk      = F0/K0;
        fkl     = log(fk);
        Kisq    = Ki.^2;
        K0sq    = K0^2; 
        
        % compute VIX
        moments.vix = (2 * er * sum(dKi.*Q./Kisq) - (fk-1)^2)/mat;
        
        % compute SVIX
        moments.svix = (2 * er * sum(dKi.*Q./K0sq) - (fk-1)^2)/mat;
        
        % compute SKEW
        eps1    = -(1+fkl-fk);
        eps2    = 2*fkl*(1-fk)+0.5*fkl^2;  % changed sign in parenth due to k/f under log
        eps3    = 3*fkl^2*(-1/3*fkl-1+fk); % changed sign in parenth due to k/f under log
        P1      = -er * sum(Q.*dKi./Kisq) + eps1;
        P2      =  er * sum(2./Kisq.*(1-log(Ki./F0)).*Q.*dKi) + eps2;
        P3      =  er * sum(3./Kisq.*(2*log(Ki./F0)-log(Ki./F0).^2).*Q.*dKi) + eps3;
        
        moments.skew   = (P3 - 3*P1.*P2 + 2*P1.^3)./(P2-P1.^2).^(3/2);      
    end
    
    if sum(strcmpi('tlm',gridparams.measnames(gridparams.meascompute)))>0
        % compute TLM
        if gridparams.tlm_limits>0.5
            goods = allputs(:,1) <= 1-gridparams.tlm_limits*sqrt(moments.mfiv/251);
        else
            goods = abs(allputs(:,5))<=gridparams.tlm_limits;
        end
        if sum(goods)>1
            Pt          = allputs(goods,3); % prices
            K           = allputs(goods,1); % strikes
            [P0,index]  = max(Pt);
            K0          = K(index);
%             [optResult  fval] = fminunc(@(x) objFunc_TailMeasure(x,P0, K0, Pt,K),[0.5 0.01]);
            [optResult,  ~, ~]   = fminsearch(@(x) objFunc_TailMeasure(x,P0, K0, Pt,K),[0.5,0.01],options);
%             fval
            moments.tlm = optResult(2)./(1-optResult(1)); % optResult(1);%
%             pause
        end
    end
    
end % if ComputeCBOEStuff

%  if isnan(moments.vix); pause; end; % for debugging -- 

end
%% tail moments 
function out=objFunc_TailMeasure(x, p0, k0, prices, strikes)
% x
xi          =   x(1);
beta        =   x(2);
checkterms  =   1+(xi/beta)*(k0-strikes);

out = sum(abs(1-p0*(checkterms.^(1-1/xi))./prices));

end
