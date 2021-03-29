function [ Erosion ] = OldMethodErosion(BeC,BeCstd,Production )
% Lupker et al. (2012) 10Be Erosion Rate Calculation

% BeC = Conc in at/g
% BeCstd = Uncertainty in at/g

% CONSTANTS
% attenuation lengths (g/mm2)
attN=160/100; 
attMs=1500/100;
attMf=4320/100;
lambda=log(2)/1.39e6;

Rau1=2.7/1000; % rock density (g/mm3)

muN=Rau1/attN;
muMs=Rau1/attMs;
muMf=Rau1/attMf;

% Get Average Production Values
rshp_Pn = reshape(Production.Pn,1,[]);
rshp_Pms = reshape(Production.Pms,1,[]);
rshp_Pmf = reshape(Production.Pmf,1,[]);
Pn =  nanmean(rshp_Pn);
Pms =  nanmean(rshp_Pms);
Pmf =  nanmean(rshp_Pmf);
PnStd = nanstd(rshp_Pn);
PmsStd = nanstd(rshp_Pms);
PmfStd = nanstd(rshp_Pmf);


% denudation rate (mm/yr)
Rate=(1./BeC).*((Pn/muN)+(Pms/muMs)+(Pmf/muMf));

%% total uncertainty estimation
n_rand=1e4; % number of random selection
means=zeros(n_rand,1); % simulated rates

% generate random concentration with truncated normal distribution
distribution = makedist('Normal');
distribution = truncate(distribution,0,inf);


% Dist A
distribution.mu = BeC;
distribution.sigma = BeCstd;
C_rand=random(distribution,n_rand,1);

% Uncertainty of neutrons: on the prod rate: 100*(0.1/3.9)=2.5% + 9% the scaling scheme
PnStd=sqrt((Pn*0.025).^2+(Pn*0.09).^2);
% Dist B
distribution.mu = Pn;
distribution.sigma = PnStd;
Pn_rand=random(distribution,n_rand,1); % random n prod rate

% Uncertainty of muon: on the prod rate: 50% + 9% the scaling scheme
PmsStd=sqrt((Pms*0.5).^2+(Pms*0.09).^2);
% DIst C
distribution.mu = Pms;
distribution.sigma = PmsStd;
Pms_rand=random(distribution,n_rand,1); % random ms prod rate
PmfStd=sqrt((Pmf*0.5).^2+(Pmf*0.09).^2);
% Dist D
distribution.mu = Pmf;
distribution.sigma = PmfStd;
Pmf_rand=random(distribution,n_rand,1); % random mf prod rate

% Compute mean values
for j=1:n_rand
    if C_rand(j)~=0
        means(j)=(1/C_rand(j))*((Pn_rand(j)/muN)+(Pms_rand(j)/muMs)+(Pmf_rand(j)/muMf));
    end
end

% generate rate pdf and find upper and lower uncertainty
rate_means=means; % simulated rates for a sample
rate_d=Rate; % rate calculated from the measured concentration and given prod rates

max_r=min(max(rate_means),25); % max simulated rate
min_r=min(rate_means); % min simulated rate)
x=min_r:0.005:max_r;
pd=fitdist(rate_means,'kernel');
y=pdf(pd,x);
mcdf=cdf(pd,x);
[~,ind]=min(abs(x-rate_d)); 
mid_cdf=mcdf(1,ind);
low_cdf=mid_cdf-0.34;
up_cdf=mid_cdf+0.34;


[~,ind]=min(abs(mcdf-low_cdf));
low_std=rate_d-x(1,ind);
[~,ind]=min(abs(mcdf-up_cdf));
up_std=x(1,ind)-rate_d;

RateStd_up=up_std;
RateStd_low=low_std;

% NEW - 05312018, just fit stuff with a lognormal dist instead of KDE
bestFitRateDist = fitdist(rate_means,'Normal');
% 

ciup = bestFitRateDist.mu + 1.96*bestFitRateDist.sigma;
cidown = bestFitRateDist.mu - 1.96*bestFitRateDist.sigma;

%%
%time scale (yr)
TimeScale=round(1./(muN*Rate));

%sediment fluxes and uncert. (Mt/yr)
Sedflux=Rau1*Production.area*Rate;
SedfluxStd=Rau1*Production.area*(RateStd_low+RateStd_up)/2;

% Add fields to struct
Erosion = struct();
Erosion.tag = Production.tag;
Erosion.area = Production.area;
Erosion.TimeScale_yr = TimeScale;
Erosion.Sedflux_mYr = Sedflux;
Erosion.SedfluxStd = SedfluxStd;
Erosion.Denudation_mmYr = Rate;
Erosion.Denudation_UpError = RateStd_up;
Erosion.Denudation_DownError = RateStd_low;
Erosion.Pn_mean = Pn;
Erosion.Pn_std = PnStd;
Erosion.Pms_mean = Pms;
Erosion.Pms_std = PmsStd;
Erosion.Pmf_mean = Pmf;
Erosion.Pmf_std = PmfStd;

% NEW
Erosion.RateSD = bestFitRateDist.sigma;
Erosion.RateMean = bestFitRateDist.mu;
Erosion.RateMeanCiUp = ciup;
Erosion.RateMeanCiDown = cidown;

end

