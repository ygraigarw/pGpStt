function C=GpPltPrm(C);

if isfield(C.RV{1},'RtrPrd')==0;
    fprintf(1,'First run GpPltTal\n. Terminating\n');
    return;
end;

% Plot parameter estimates nicely
clf;
tTtl={'QR';'Pss';'GP shape';'GP scale';sprintf('%g-yr RV',C.RV{1}.RtrPrd)};
for iN=1:C.nNep;
    for j=1:5;
        switch j
            case 1; tDat=C.QR{iN}.Prm(C.QR{iN}.n2Plt+1:end,1); tTxt=C.QR{iN}.PrmNms{1};
            case 2; tDat=C.Pss{iN}.Prm(C.Pss{iN}.n2Plt+1:end,1); tTxt=C.Pss{iN}.PrmNms{1};
            case 3; tDat=C.GP{iN}.Prm(C.GP{iN}.n2Plt+1:end,1); tTxt=C.GP{iN}.PrmNms{1};
            case 4; tDat=C.GP{iN}.Prm(C.GP{iN}.n2Plt+1:end,2); tTxt=C.GP{iN}.PrmNms{2};
            case 5; tDat=C.RV{iN}.RV; tTxt='RV';
        end;
        subplot(2,3,j); hold on;
        pDns(tDat,10,pClr(iN));
        if iN==C.nNep;
            title(tTtl{j});
            xlabel(tTxt,'interpreter','LaTeX');
            pAxsLmt; pDflHug;
        end;
    end;
end;
tClr={'r';'m';'g';'c';'b';'k';'gr';'br'};
t='NEPs: ';
for iN=1:C.nNep;
    t=sprintf('%s %g (%s)',t,C.Nep(iN),tClr{iN});
end;
pDatStm(t);
pGI('DgnPltPrm',2);

return;

