%% Simple code to fit QR, Pss and GP models, all the parameters assumed to be fixed in time
% Mdl(xi, sgm, rho, psi) with xi=xi_0, sgm=sgm_0, ...
% Will run as is for toy data to check 
% Need to input data as structure (see occurrences of USER INPUT) below
%
% QR = non-stationary Quantile Regression threshold
% Pss = non-stationary PoiSSon could model for threshold exceedances
% GP = Generalised Pareto model for size of threshold exceedances
%
% This is a simplified version of pGpNonStt

%% Output for further plotting / investigation etc
% Output from the analysis is saved in a structure C
% A typical structure C (when 8 different threshold non-exceedance probabilities are used) is
%
%% C
% 
%       Nep: [8×1 double] Non-exceedance probabilities for thresholds 
%      nNep: 8            Number of NEPs
%        QR: {[1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]}  QR model details (inc posterior sample)
%       Pss: {[1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]}  Pss model details (inc posterior sample)
%        GP: {[1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]}  GP model details (inc posterior sample)
%        RV: {[1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]}  RV details (inc posterior sample)
%     RVCmp: [1×1 struct]                                                                                                      RV comparison summary
%    PrmSmm: [1x1 struct]                                                                                                      Assessment of slope parameter changes
%
%% C.QR{q}, C.Pss{q}, C.GP{q} for q=1,2,..., nNep are structures like
% 
%        Lkl: 'QR'
%       nPrm: 1
%     PrmNms: {1×1 cell}
%        Nep: 0.6000
%       nItr: 10000
%      n2Plt: 5000
%     NgtStr: 0.1000
%     AdpItr: 1000
%     AdpBet: 0.0500
%     PrmStr: [2×1 double]
%     AccRat: [10000×1 double]
%        Prm: [10000×2 double]
%        Nll: [10000×1 double]
%     PrmUse: [1×1 double]
%
%% Key output for further plotting etc are
%
% C.QR{q}.Prm   nItr x 1 values of psi0 from MCMC (in general it is safe to use the last 9000; first 1000 might involve "burn-in")
% C.Pss{q}.Prm  nItr x 1 values of rho0 from MCMC (in general it is safe to use the last 9000; first 1000 might involve "burn-in")
% C.GP{q}.Prm   nItr x 2 values of xi0, sigma0 from MCMC (in general it is safe to use the last 9000; first 1000 might involve "burn-in")
% C.RV{q}.RV    nRls x 1 values of return values generated using C.QR.Prm, C.Pss.Prm and C.GP.Prm
% C.RVCmp.Qnt   nNep x 3 Quantlies (2.5%, 50% and 97.5%) at "time end" and quantlies (2.5%, 50% and 97.5%) at "time end" for each NEP

%% Set up
clc; clear; clf;
VrbNms={'$\xi$';'$\sigma$';'$\psi$'};

%% Simulate a sample of data
if 1; %for testing

    % True parameters P0=[xi0;sgm0;rho0;psi0;] of linear regression 
    
    %% *** USER INPUT *** Pick the type of simulated data
    X.Prm0=[-0.3; 2; 20; 2;];
    
    % Time variable
    % NB A COMMON time value used for observations in the same year
    X.nYr=85;
    tYr=(1:X.nYr)';  % Time in years

    % True parameter estimates per year
    X.XSM0=[ones(X.nYr,1)*X.Prm0(1) ones(X.nYr,1)*X.Prm0(2) ones(X.nYr,1)*X.Prm0(3) ones(X.nYr,1)*X.Prm0(4)];

    % Number of occurrences per annum
    tOcc=poissrnd(X.XSM0(:,3));
    
    % Generate data from GP
    k=0;
    X.nT=sum(tOcc);
    X.Tim=nan(X.nT,1);
    X.Dat=nan(X.nT,1);
    for iY=1:X.nYr;
        for iO=1:tOcc(iY);
            k=k+1;
            X.Tim(k)=tYr(iY)/X.nYr;
            X.Dat(k)=gprnd(X.XSM0(iY,1),X.XSM0(iY,2),X.XSM0(iY,4));
        end;
    end;
    
    X, % See the structure
    
    subplot(2,2,1); plot(tYr,tOcc,'ko');
    subplot(2,2,2); plot(X.Tim*X.nYr,X.Dat,'ko');
    
end;

%% ***USER INPUT*** Read in your data here
if 0;
    %X.nYr ; %  1   x 1 number of years
    %X.nT  ; %  1   x 1 number of occurrences
    %X.Tim ; % nT   x 1 years on [0,1] (so that floor((X.Tim*X.nYr)+1) gives the year number
    %X.Dat ; % nT   x 1 data
    %     load('G:\UoM Climate Change\ssp245.mat');
    Fld=Field;
    %     X.Dat=data(:,3);
    %     t1=data(:,1);
    X.Dat=POT.(Fld)(:,2);
    t1=POT.(Fld)(:,1);
    t2=(floor(t1(1)):floor(t1(end)))';
    X.Tim = (floor(t1)-floor(t1(1))+1)/(floor(t1(end))-floor(t1(1))+1);
%     X.Tim=(t2-t2(1))/range(t2);
%     X.Tim=1:numel(t1)/numel(t1);
    X.nT=size(t1,1);
    X.nYr = numel(unique(floor(t1)));
end;

%% ***USER INPUT*** Read in your data here - MADAGASCAR ANALYSIS
if 0;
     load Madagascar.mat;
    
     nT=size(yrs,1);
     Tim=yrs;
     Dat=HsPOTall;
     
     % Kevin's code to create structure X, adapted from above
     t1=yrs;
     t2=(floor(t1(1)):floor(t1(end)))';
     X.Tim = (floor(t1)-floor(t1(1))+1)/(floor(t1(end))-floor(t1(1))+1);
     X.nT=size(t1,1);
     X.nYr = numel(unique(floor(t1)));
     X.Dat=Dat;
     
     plot(X.Tim,X.Dat,'k.');
     
end;

%% ***USER INPUT*** Specify NEPs to consider
if 1;
    %C.Nep=(0.6:0.05:0.95)'; % (0.6:0.05:0.9)' is a good range; but maybe you want to use (0.7:0.1:0.9)' to get going
    %C.Nep=(0.7:0.1:0.9)';
    C.Nep=[0.9;(0.95:0.01:0.99)';0.995];
    %
    C.nNep=size(C.Nep,1);
end;

%% Estimate extreme value threshold (linear Quantile Regression)
if 1;
    
    for iN=1:C.nNep
        
        C.QR{iN}.Lkl='QR';       % Likelihood
        C.QR{iN}.nPrm=1;         % Number of parameters
        C.QR{iN}.PrmNms={'$\psi_0$';}; % Names for parameters
        C.QR{iN}.Nep=C.Nep(iN);  % NEP

        C.QR{iN}.nItr=10000;     % Number of MCMC iterations - 1e4 minipsim when used in anger
        C.QR{iN}.n2Plt=5000;     % Number of iterations from end of chain to "beleive"
        
        C.QR{iN}.NgtStr=0.1;     % Candidate random walk standard deviation - don't change
        C.QR{iN}.AdpItr=1000;    % Number of warm up iterations - don't change
        C.QR{iN}.AdpBet=0.05;    % Adaptive MC - don't change        C.Nep=X.Nep(iN);
        
        C.QR{iN}.PrmStr=[quantile(X.Dat,C.Nep(iN))]; % Constant starting solution for quantile regression
        
        C.QR{iN}=McmcStt(X,C.QR{iN});   % Run MCMC algorithm
        
        C.QR{iN}.PrmUse=mean(C.QR{iN}.Prm(C.QR{iN}.nItr-C.QR{iN}.n2Plt+1:C.QR{iN}.nItr,:))'; % Use posterior mean for subsequent inference
        
        tStr=sprintf('Mdl%s-Nep%g',C.QR{iN}.Lkl,C.QR{iN}.Nep); pDatStm(tStr); pGI(tStr,2); % Save plot  
        tFil=sprintf('MCMC'); save(tFil,'C'); % Save whole chain

    end;
    
end;

%% Estimate rate of threshold exceedance per annum (linear Poisson Process)
if 1;
        
    for iN=1:C.nNep
        
        C.Pss{iN}.Lkl='Pss';      % Likelihood
        C.Pss{iN}.Nep=C.Nep(iN);  % NEP
        C.Pss{iN}.PrmNms={'$\rho_0$';}; % Names for parameters
        C.Pss{iN}.nPrm=1;         % Number of parameters
        
        C.Pss{iN}.nItr=10000;     % Number of MCMC iterations - 1e4 minipsim when used in anger
        C.Pss{iN}.n2Plt=5000;     % Number of iterations from end of chain to "beleive"
        
        C.Pss{iN}.NgtStr=0.1;     % Candidate random walk standard deviation - don't change
        C.Pss{iN}.AdpItr=1000;    % Number of warm up iterations - don't change
        C.Pss{iN}.AdpBet=0.05;    % Adaptive MC - don't change        C.Nep=X.Nep(iN);
        
        % Estimate Poisson count for threshold exceedances
        t1=(X.Dat-ones(X.nT,1)*C.QR{iN}.PrmUse(1))>0; % Threshold exceedances
        for iY=1:X.nYr;
            t2=floor(X.Tim*X.nYr)>=(iY-1) & floor(X.Tim*X.nYr)<iY; % Particular year
            C.Pss{iN}.Cnt(iY,:)=sum(t1(t2==1));
            C.Pss{iN}.CntTim(iY,:)=(iY-0.5)/X.nYr; % Take middle of year
        end;
        
        % Constant starting solution from Poisson fit
        C.Pss{iN}.PrmStr=[poissfit(C.Pss{iN}.Cnt)]; 
        
        C.Pss{iN}=McmcStt(X,C.Pss{iN});   % Run MCMC algorithm
        
        tStr=sprintf('Mdl%s-Nep%g',C.Pss{iN}.Lkl,C.Pss{iN}.Nep); pDatStm(tStr); pGI(tStr,2); % Save plot
        tFil=sprintf('MCMC'); save(tFil,'C'); % Save whole chain

    end;
    
end;

%% Estimate size of threshold exceedance per annum (linear Generalised Pareto)
if 1;
            
    for iN=1:C.nNep
        
        C.GP{iN}.Lkl='GP'; 
        C.GP{iN}.Nep=C.Nep(iN);  % NEP
        C.GP{iN}.PrmNms={'$\xi_0$';'$\sigma_0$';}; % Names for parameters
        C.GP{iN}.nPrm=2;     % Number of parameters
        
        C.GP{iN}.nItr=10000;     % Number of MCMC iterations - 1e4 minipsim when used in anger
        C.GP{iN}.n2Plt=5000;     % Number of iterations from end of chain to "beleive"
        
        C.GP{iN}.NgtStr=0.1;     % Candidate random walk standard deviation - don't change
        C.GP{iN}.AdpItr=1000;    % Number of warm up iterations - don't change
        C.GP{iN}.AdpBet=0.05;    % Adaptive MC - don't change        C.Nep=X.Nep(iN);
        
        % Isolate threshold exceedances and times of occurrence
        t1=X.Dat-(ones(X.nT,1)*C.QR{iN}.PrmUse(1)); % Threshold exceedances
        C.GP{iN}.Exc=t1(t1>0);
        C.GP{iN}.ExcTim=X.Tim(t1>0);
        
        % Constant starting solution from GP fit
        t=gpfit(C.GP{iN}.Exc); 
        if t(1)<-0.5; t(1)=-0.4; end;
        if t(1)>0.5; t(1)=0.4; end;
        C.GP{iN}.PrmStr=[t(1);t(2)]; 
        
        C.GP{iN}=McmcStt(X,C.GP{iN});   % Run MCMC algorithm
        
        tStr=sprintf('Mdl%s-Nep%g',C.GP{iN}.Lkl,C.GP{iN}.Nep); pDatStm(tStr); pGI(tStr,2); % Save plot  
        tFil=sprintf('MCMC'); save(tFil,'C'); % Save whole chain

    end;
    
end;

%% Plot tails per threshold, RV with threshold, and parameter estimates nicely
if 1;
    
    RtrPrd=100; % return period to use in years
    nRls=1000; % number of realisations of tails to generate
    C=GpPltTal(C,X,RtrPrd,nRls); % plot of tails per threshold, and RV quantiles with threshold
    C=GpPltPrm(C); % nice plot of QR-Pss-GP parameters for all thresholds

end;

%% Update output file
if 1;
    
    tFil=sprintf('MCMC'); save(tFil,'C'); % Save whole chain

end;