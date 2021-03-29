function S_Z= oneD_TVD_eroder(S_Z, n, dt, A, d, r, dx, dd, rr, dx_centered, U,cfls)
% Explicit solution for river incision, calucalted by solving the Stream Power
% Law.
%
% Syntax
%
%       Z = oneD_TVD_eroder(S_Z, n, dt, A, d, r, dx, dd, rr, dx_centered, U)
%
% Description
%
%       Explicit solution for river incision, calucalted by solving the
%       Stream Power Law. Internally, time steps are adapted as to honor
%       the cfl criterion.
%
% Input
%
%       p         parameter values
%       Z         digital elevation model (matrix)
%       dt        time step
%       A         velocity field
%       d         donor, derived from the flow
%                 direction FD.ix where ix is an edge attribute and
%                 represetns topologically sorted nodes (givers)
%       r         receiver, derivedfrom the flow direction FD.ixc where
%                 ixc is an edge attribute and represetns topologically
%                 sorted nodes (receivers) dx_ik horizontal distance
%       dx        distance between giver and receiver
%       dd        donor of donor 
%       rr        Receiver of receiver 
%       dx_centered distance between giver of giver and receiver
%       upl       Uplift rate (m/yr)
%       cfls      Stability number
%
% Output
%
%       Z       digital elevation model with adapted river elevation
%
% Example
%
%
% See also:
%
% Authors: Benjamin Campforts (benjamin.campforts[at]ees.kuleuven.be)
%          Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%
% Date: 28. Januari, 2015



a_ori=-A(d);
a_m = min(0,a_ori);
a_p = max(0,a_ori);
a_vect=a_ori;

nrtsteps = ceil(dt/(cfls*min(dx)/max(abs(a_vect))));
dte = dt/nrtsteps;
time=dt;

while time>0;
    %% TVD
    % Nonlinear solution, calculated by integrating the explicit FDM
    % solution to the incision speed
    if n~=1
        s1=((max(S_Z(d)-S_Z(r),0))./dx); %dx_ik? 
        %Very low slopes are leading through very small timesteps so they are set to 0.
        s1(s1<1e-4)=0;
        exp_f=s1.^(n-1);
        exp_f(isinf(exp_f))=1;
        a_vect(exp_f~=0)=a_ori(exp_f~=0).*exp_f(exp_f~=0);
        a_m = min(0,a_vect);
        a_p = max(0,a_vect);
        
        % Check timestep
        dt_calc = cfls*min(dx)/max(abs(a_vect));
        if dt_calc<dte
            disp(['For stability, dte is set to: ' num2str(dt_calc)]);
            dte=dt_calc;
        end
    end
    %Time
    time=time-dte;
    if time<0
        dte=dte+time;
        time=0;
    end
       
   
        el_c=S_Z(d);
        el_d=S_Z(r);
        kk_nb=rr;
            kk_nb(isnan(rr))=1;
            el_d2=S_Z(kk_nb);    
            el_d2(isnan(rr))=nan;%el_d(isnan(kk));
            iiNb=dd;
            iiNb(isnan(dd))=1;
            el_up=S_Z(iiNb);
            el_up(isnan(dd))=nan;
        r_TVD=(el_d2-el_d)./(el_d-el_c);
        r_TVD((el_d-el_c)==0)=1;
        
        % Define Flux Limiter function
        %VANLEER
        phi = (r_TVD + abs(r_TVD))./(1 + abs(r_TVD));
        
        l_TVD=(el_d-el_c)./(el_c-el_up);
        l_TVD((el_c-el_up)==0)=1;
        
        % Define Flux Limiter function
        %VANLEER
        phi_l = (l_TVD + abs(l_TVD))./(1 + abs(l_TVD));        
        
        % Compute fluxes for TVD
        F_rl = a_p.*el_c + a_m.*el_d;
        F_rh = (1/2)*a_vect.*(el_c+el_d) - (1/2)*(a_vect.^2).*...
            (dte./dx_centered).*(el_d-el_c);
        F_ll = a_p.*el_up + a_m.*el_c;
        F_lh= (1/2)*a_vect.*(el_up+el_c) - (1/2)*(a_vect.^2).*...
            (dte./dx_centered).*(el_c-el_up);
             
        F_right = F_rl + phi.*(F_rh-F_rl);
        F_left = F_ll+ phi_l.*( F_lh- F_ll);
        TVD_next= S_Z(d)-(dte*(F_right-F_left)./dx_centered);
        TVD_next(TVD_next<0)=0;
        
        % In case of nan value, replace by explicit solution        
        exp_next=S_Z(d)-dte*A(d).*((max(S_Z(d)-S_Z(r),0))./dx);        
        TVD_next(isnan(TVD_next))=exp_next(isnan(TVD_next));        
        S_Z(d) = TVD_next;        
        S_Z=S_Z+dte.*U;
      
end
end
