function outmat = OM_Parallel4DataHandling(optdata0, zerocd, gridparams)
% version 2013-01-01
% Please write to vilkov@vilkov.net with questions/comments
% 2014 Grigory Vilkov

entrycount = 0;
outmat     = NaN(252*gridparams.matnum,length(gridparams.measnames)+3);

if size(optdata0,1)>0
    % go thru each date for a given secid
    dates = unique(optdata0(:,2));
    for i=1:length(dates)
%         disp(i)
        optdata1 = optdata0(optdata0(:,2)==dates(i),:);
        
        % interpolate zerocd to get the rate
        zerocurr = zerocd(zerocd(:,1)==dates(i),2:3); % select zerocd for a given date
        if size(zerocurr,1)>1
            currspline = pchip(zerocurr(:,1),zerocurr(:,2)); % estimate current spline
        elseif size(zerocurr,1)==1
            currspline = pchip([1; 100],zerocurr(:,2).*ones(2,1));
        else
            currspline = pchip([1; 100],[1; 1]);
        end
        % check the available maturities for the given secid/date combo
        days = unique(optdata1(:,3));
        
        % go through each maturity separately
        for q=1:length(days)
            % selected all options with a maturity needed, now we can proceed computing the moments, etc.
            currsel     = optdata1(:,3)==days(q) & ~isnan(sum(optdata1,2));
            optdata2    = optdata1(currsel,:);
            
            % if there are more than 3 options, proceed to computing
            if sum(currsel)>3
                % evaluate rate
                rate  = ppval(currspline,days(q));
                if rate>=1; rate=rate/100; end; % just in case we have rates in pct
                % ###################
                % COMPUTE IMPLIED MOMENTS                               
                [moments, nopt]=OM_Compute_Implied_Moments(optdata2, rate, gridparams);
                % and write them to the output matrix
                entrycount               = entrycount+1;
                for z=1:length(gridparams.measnames)
                    outmat(entrycount,z)    = moments.(gridparams.measnames{z});
                end
                outmat(entrycount,end-2:end)    = [dates(i) days(q) nopt];
            end
        end
        
    end
end


end