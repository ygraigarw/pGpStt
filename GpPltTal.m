function C=GpPltTal(C,X,RtrPrd,nRls);

[nr,nc]=pNmbSubPlt(C.nNep+1);

%% Calculate tails per threshold choice, RV variability per threshold
for iN=1:C.nNep;
    
    C.RV{iN}.RtrPrd=RtrPrd; % Return period of interest
    C.RV{iN}.nRls=nRls;  % Number of realisations to use (1000 is good)
    
    % Parameter estimates at start
    t1=randi(C.QR{iN}.nItr-C.QR{iN}.n2Plt,C.RV{iN}.nRls,1)+C.QR{iN}.n2Plt;
    t2=randi(C.Pss{iN}.nItr-C.Pss{iN}.n2Plt,C.RV{iN}.nRls,1)+C.Pss{iN}.n2Plt;
    t3=randi(C.GP{iN}.nItr-C.GP{iN}.n2Plt,C.RV{iN}.nRls,1)+C.GP{iN}.n2Plt;
    %For testing: tPsi=C.QR{iN}.PrmUse(1); %This uses the MEAN threshold (and so does not propagate threshold uncertainty)
    tPsi=C.QR{iN}.Prm(t1,1);%This uses a random threshold (and so "slightly" overestimates the effect of threshold uncertainty)
    tRat=C.Pss{iN}.Prm(t2,1); % Rate is per annum, modelled at midpoint of year
    tXi=C.GP{iN}.Prm(t3,1);
    tSgm=C.GP{iN}.Prm(t3,2);
    
    C.RV{iN}.RV(:,1)=(tSgm./tXi).*( (C.RV{iN}.RtrPrd.*tRat).^tXi - 1) + tPsi; % Return value
    
    % Summary statistics of differences
    C.RVCmp.Qnt(iN,:)=[quantile(C.RV{iN}.RV(:,1),[0.025 0.5 0.975])]; % Quantiles
    
    % Plot
    subplot(nr,nc,iN); hold on;
    
    % Generate data from QR-Pss-GP model and plot tails
    nTal=100;
    tnYr=size(C.Pss{1}.Cnt,1);
    for iR=1:nTal;
        iI=randi(size(tRat,1));
        tn=round(tRat(iI)*tnYr);
        tDat=gprnd(tXi(iI)*ones(tn,1),tSgm(iI),tPsi(iI));
        plot(sort(tDat),log10(1-(((1:tn)'-0.5)/tn))+log10(1-C.QR{iN}.Nep),'k.')
    end;
    te=size(C.GP{iN}.Exc,1);
    plot(sort(C.GP{iN}.Exc)+tPsi(iI),log10(1-(((1:te)'-0.5)/te))+log10(1-C.QR{iN}.Nep),'r*')
    pAxsLmt; pDflHug;
    title(sprintf('Tail: NEP=%g',C.Nep(iN)));
    
end;

% Plot return value quantiles
subplot(nr,nc,nr*nc);hold on;
plot(C.Nep,C.RVCmp.Qnt(:,1),'ko-');
plot(C.Nep,C.RVCmp.Qnt(:,2),'ko-','linewidth',3);
plot(C.Nep,C.RVCmp.Qnt(:,3),'ko-');

if isfield(X,'Prm0')==1; % True values are known
    tPsi=X.Prm0(4);
    tRat=X.Prm0(3);
    tSgm=X.Prm0(2);
    tXi=X.Prm0(1);
    tRVTru=(tSgm./tXi).*( (C.RV{1}.RtrPrd.*tRat).^tXi - 1) + tPsi;
    plot(C.Nep,ones(C.nNep,1)*tRVTru,'r--');
end;
pAxsLmt; pDflHug;
title 'Quantiles of RV [True=red]'; xlabel 'Threshold NEP';

pGI('DgnPltTal',2);

return;
