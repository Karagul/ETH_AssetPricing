
% ####################################################
% COMPUTE IMPLIED MOMENTS: MAIN ROUTINE
% USING OWN ROUTINE WITH INTERP: MFIV(variance)/MFIS/MFIK/CXI(variance) (corridor MFIV with interpolation)
% USING CBOE ROUTINE:            VIX(as variance)/ SKEW
% version July 01, 2014
% Please write to vilkov@vilkov.net with questions/comments
% 2014 Grigory Vilkov


% The input files: 
% 1. OPTIONS: opprcd_4mf2000INDX_MAIN contains RAW (opprcd) data for year 2000
% for several indices as follows: 
% secid ticker
% 101507 MID
% 102442 SML
% 102456 DJX
% 108105 SPX
% 109764 OEX
% 2. ZEROCD: zerocd - zero CD rates
% 3. ID and DATES: id_dt

% ## OPTIONS: INPUT FILE
% 1 b.secid,
% 2 b.date,
% 3 b.exdate-b.date as days,
% 4 b.last_date,
% 5 b.strike_price/abs(b1.close) as mnes, (beware of 1000x strike in OptionMetrics)
% 6 b.best_bid/abs(b1.close),       
% 7 b_best_offer/abs(b1.close),
% 8 b.volume,
% 9 b.open_interest,
% 10 b.impl_volatility,
% 11 b.delta (from -1 to 1)

% ## FILTERS FOR OPTIONS DATA - use the switch below: gridparams.filter  
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
%         prc/spot<0.1%


% ### OUTPUT: structure IM that contains a number of vectors with the mf measures
% the last three fields contain the information of the date on which a
% given measure was computed, number of options used for the routine with
% interpolation, and days to maturity. 
% IM structure can be then allocated to the matrices as needed.
% IM = 
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
%     'iv'            %12 IV
%       dt  
%     nopt  
%     days  


clear all
tic


% ###  DEFINE CONTROLS ####

% as input use annual files 
YrFrom                  = 2000
YrTo                    = 2000

NumWorkers = 0; % set to 0 if no parallel, set to # of workers for parallel


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
    'iv'   }        %12 ATM IV                       ()


% grid/data parameters for implied moments
gridparams.k            = 2;   % limits are from 1/(1+k) to 1+k
gridparams.m            = 500; % use 2*m+1 points, i.e. -m:m
gridparams.interp_flag  = 1;
gridparams.matnum       = 10; % number of maturities for each secid
gridparams.measnames    = OIMeasureNames(:,1);
gridparams.meascompute  = [1 2 3 4 5 7 8 10 12]; % 
gridparams.filter       = 2; % filter for 'our' routine; the ther routine uses CBOE rules/ NOT USED -- ONLY IN FILE NAMES
gridparams.mnes_limits  = [0.5 1.5];        % filter out the options beyond these mnes limits, by default [0.7 1.3]

gridparams.cx_limits    = [0.03 0.97];
gridparams.tlm_limits   = 1.65;             % if <0.5, then delta or otherwise # sigmas limit for the tail measure
gridparams.matlim       = [7 300];          % maturities for each secid
gridparams.surface      = false;            % set to true if bids=ask 

MatNum   = gridparams.matnum;
dataflds = {OIMeasureNames{:,1},  'dtt', 'days', 'nopt'};


%%
% ###  DEFINE FILE NAMES - LOAD SOME DATA ####
% sample
SMPL     = 'INDX_MAIN';

% folders
main_dir1 = '';
mat_dir2  = '';

%  mat files
fn1.fn_secid_dt              = 'id_dt';
fn1.fn_zerocd                = 'zerocd';
fn1.fn_data                  = 'opprcd_4mf';
fn1.fn_implied_mommat        = 'mf_moments';

flds  = fields(fn1);
for z=1:length(flds)
    fn1.(flds{z}) = [main_dir1 mat_dir2 fn1.(flds{z})];
end


% #### parallel
% CurrWorkers=parpool('size');
% if  CurrWorkers~=NumWorkers && CurrWorkers>0 && NumWorkers>0
%     parpool close force
%     parpool('open', 'local', NumWorkers)
% elseif NumWorkers>0
%     parpool('open', 'local', NumWorkers)
% end

%%

% #########################################################################
% ###### COMPUTE IMPLIED MOMENTS
% load the necessary files
load(fn1.fn_secid_dt)
load(fn1.fn_zerocd);

dtall = dtall(:,1);
% ### define and apply filters


for yr=YrFrom:YrTo
    
    % load the data
    fn_datamat = [fn1.fn_data num2str(yr) '' SMPL '.mat'];    
    load(fn_datamat);
    
        IM.mfiv = [];
        IM.smfiv= [];
        IM.mfis = [];
        IM.mfik = [];
        IM.vix  = [];
        IM.svix = [];
        IM.skew = [];
        IM.cx   = [];
        IM.cxi  = [];
        IM.tlm  = [];
        IM.tlmi = [];
        IM.iv   = [];
        IM.dt   = [];
        IM.nopt = [];
        IM.days = [];
        
        mfiv = NaN(252*MatNum,length(secid));
        smfiv= NaN(252*MatNum,length(secid));
        mfis = NaN(252*MatNum,length(secid));
        mfik = NaN(252*MatNum,length(secid));
        vix  = NaN(252*MatNum,length(secid));
        svix = NaN(252*MatNum,length(secid));
        skew = NaN(252*MatNum,length(secid));
        cx   = NaN(252*MatNum,length(secid));
        cxi  = NaN(252*MatNum,length(secid));
        tlm  = NaN(252*MatNum,length(secid));
        tlmi = NaN(252*MatNum,length(secid));
        iv   = NaN(252*MatNum,length(secid));
        dtt  = NaN(252*MatNum,length(secid));
        days = NaN(252*MatNum,length(secid));
        nopt = NaN(252*MatNum,length(secid));
        

               
        % remove shorter than LIM and loner than LIM days
        bads            = optdata(:,3)<gridparams.matlim(1) | optdata(:,3)>gridparams.matlim(2);
        optdata(bads,:) = [];
        
        % go thru each secid, either as parallel or serial
        parfor j=1:length(secid)
            
            currsel  = optdata(:,1)==secid(j);
            optdata0 = optdata(currsel,:);
            
            outmat = OM_Parallel4DataHandling(optdata0, zerocd, gridparams);
            
            mfiv(:,j) = outmat(:,1);
            smfiv(:,j)= outmat(:,2);
            mfis(:,j) = outmat(:,3);
            mfik(:,j) = outmat(:,4);
            cxi(:,j)  = outmat(:,5);
            tlmi(:,j) = outmat(:,6);
            vix(:,j)  = outmat(:,7);
            svix(:,j) = outmat(:,8);
            skew(:,j) = outmat(:,9);
            cx(:,j)   = outmat(:,10);
            tlm(:,j)  = outmat(:,11);
            iv(:,j)   = outmat(:,12);
            dtt(:,j)  = outmat(:,13);
            days(:,j) = outmat(:,14);
            nopt(:,j) = outmat(:,15);
        end
        
        
        %             data flds
        
        bads = (nansum([mfiv vix],2))==0;
        mfiv(bads,:) = []; IM.mfiv = [IM.mfiv; mfiv];
        smfiv(bads,:)= []; IM.smfiv= [IM.smfiv; smfiv];
        mfis(bads,:) = []; IM.mfis = [IM.mfis; mfis];
        mfik(bads,:) = []; IM.mfik = [IM.mfik; mfik];
        vix(bads,:)  = []; IM.vix  = [IM.vix;   vix];
        svix(bads,:) = []; IM.svix = [IM.svix;   svix];
        skew(bads,:) = []; IM.skew = [IM.skew; skew];
        cx(bads,:)   = []; IM.cx   = [IM.cx;     cx];
        cxi(bads,:)  = []; IM.cxi  = [IM.cxi;   cxi];
        tlm(bads,:)  = []; IM.tlm  = [IM.tlm;   tlm];
        tlmi(bads,:) = []; IM.tlmi = [IM.tlmi; tlmi];
        iv(bads,:)   = []; IM.iv   = [IM.iv;     iv];
        dtt(bads,:)  = []; IM.dt   = [IM.dt  ;  dtt];
        days(bads,:) = []; IM.days = [IM.days; days];
        nopt(bads,:) = []; IM.nopt = [IM.nopt; nopt];
    
    
    disp(['Implied Moments computed '  '/ Year: ' num2str(yr) ' Time: ' num2str(toc)])
    
    clear optdata*;
    
    fn_implied_mommatX   =[fn1.fn_implied_mommat '_year' num2str(yr) '_filter' num2str(gridparams.filter)];
    save(fn_implied_mommatX,'IM')
    
end








