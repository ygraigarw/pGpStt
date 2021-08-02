function Nll=EstNllStt(Prm, Y, C);
%function Nll=EstNllStt(Prm, Y, C);
%
%Negative log likelihood for QR, Pss, GP and GEV, All parameters are constant in time T

X=Y.Dat;
T=Y.Tim;

Nll=NaN;

switch C.Lkl;
    
    case 'QR';
        
        %% Linear forms for parameters
        Psi=Prm(1);
        
        %% Calculate sample NLL
        t0=X-Psi;
        t(t0>0)=t0(t0>0).*C.Nep;
        t(t0<=0)=-t0(t0<=0)*(1-C.Nep);
        Nll=sum(t);
        
    case 'Pss';
        % Note that we work with annual counts here
        
        %% Linear forms for parameters
        Rho=ones(size(C.Cnt,1),1)*Prm(1);
        
        %% Reject negative rate
        if sum(Rho<0)>0;
            return;
        end;
        
        %% Calculate sample NLL
        d=1; % Annual rate
        Nll=d*sum(Rho)-sum(C.Cnt.*log(Rho));
        
    case 'GP';
        
        %% Linear forms for parameters
        Xi=Prm(1);
        Sgm=Prm(2);
        
        %% Reject non-positive sigma
        if min(Sgm)<0;
            return;
        end;
        
        %% Reject xi out of sensible range
        if min(Xi)<=-0.5 || max(Xi)>0.5;
            return;
        end;
        
        %% Reject exceedances of upper end point
        t0=(1 + Xi.*C.Exc./Sgm);
        if min(t0)<0;
            return;
        end;
        
        %% Calculate NLL
        t0=(1 + Xi.*C.Exc./Sgm);
        Nll = sum(log(Sgm)+(1./Xi+1).*log(t0));
        
end;

return;