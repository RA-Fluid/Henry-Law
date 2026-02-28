%--------------------------------------------------------------------------
% This code is developed to calculate the
%
% - Henry's volatility constants, Hv_px,
% - Henry's solubility constants, Hs_cp
% - molar fractions,
% - molarity,
% - mass fractions,
% - partial molar volumes, and
%- variation in the density of water by concentration of each gas component
% 
% based on the Henry's laws & some other sources. Please refer to the
% report for the details
%
% Developed by Reza Azadi, PhD at the University of Sydney (2026)
%
%--------------------------------------------------------------------------

%% initialise

clc
close all
clear

% Add dependencies to the MATLAB path
addpath(genpath('functions'));

%% Main inputs
% [N2, O2, Ar]

% volume fractions of nitrogen (N2), oxygen (O2), Argon (Ar), and carbon dioxide (CO2) in air
volf = [0.78103, 0.20940, 0.00917, 0.00040] ; % []

% molar mass of nitrogen (N2) and oxygen (O2)
molm = [28.01, 32.00, 39.948, 44.01] ; % [gr mol^-1]
% molar mass of water (H2O)
molm_f = 18.01528 ; % [gr mol^-1]

% gas constant
R = 8.3145 ; % [m^3 Pa K^-1 mol^-1]

% atm pressure
p_atm = 101.325 ; % [kPa]

%% Define the temperature and pressure ranges

% total pressure
P_abs = ( (0:100:1000) + p_atm )' ; % [kPa]
Temp = (20 : 2 : 30)' ; % [Celsius]

NT = numel(Temp) ;
NP = numel(P_abs) ;

% water saturation pressure at the selected temperatures
Pw_sat = zeros(1, NT) ;
Pt = zeros(NP, NT) ;
for i = 1 : NT
    Pw_sat(i) = XSteam('psat_T', Temp(i)) * 100 ; % [kPa] water

    for j = 1 : NP
        Pt(j,i) = P_abs(j) - Pw_sat(i) ; % total pressure corrected by the water sta pressure at the desired temp [kPa]
    end
end

%% Calculate the density of the solvent (water here)

densf = zeros(NP, NT) ;
for i = 1 : NP
    for j = 1 : NT
    densf(i,j) = XSteam('rho_pT', Pt(i)*0.01, Temp(j)) ; % [kg m^-3] water
    end
end

%% Calculate the fractions
% partial pressures of nitrogen (N2), oxygen (O2), Argon (Ar) and CO2 in air

Pg_N2 = volf(1) .* Pt ; % [kPa]
Pg_O2 = volf(2) .* Pt ; % [kPa]
Pg_Ar = volf(3) .* Pt ; % [kPa]
Pg_CO2 = volf(4) .* Pt ; % [kPa]

%% Calculate Henry's volatility constant, Hv_px

Hv_px_Pa = zeros(4, NT) ;
for i = 1: numel(Temp)
    Hv_px_Pa(1,i) = RA_HenryLaw('nitrogen', Temp(i)+273.15) ; % [Pa]
    Hv_px_Pa(2,i) = RA_HenryLaw('oxygen', Temp(i)+273.15) ; % [Pa]
    Hv_px_Pa(3,i) = RA_HenryLaw('Argon', Temp(i)+273.15) ; % [Pa]
    Hv_px_Pa(4,i) = RA_HenryLaw('carbon dioxide', Temp(i)+273.15) ; % [Pa]
end

Hv_px_atm  = Hv_px_Pa ./ 101325 ; % [atm]

%% Calculate Henry's solubility constant, Hs_cp

Hs_cp_N2 = zeros(size(densf)) ;
Hs_cp_O2 = zeros(size(densf)) ;
Hs_cp_Ar = zeros(size(densf)) ;
Hs_cp_CO2 = zeros(size(densf)) ;

for i = 1 : NT
    Hs_cp_N2(:,i) = (densf(:,i) ./ (molm_f*1e-3)) / Hv_px_Pa(1,i) ; % [mol m^-3 Pa^-1]
    Hs_cp_O2(:,i) = (densf(:,i) ./ (molm_f*1e-3)) / Hv_px_Pa(2,i) ; % [mol m^-3 Pa^-1]
    Hs_cp_Ar(:,i) = (densf(:,i) ./ (molm_f*1e-3)) / Hv_px_Pa(3,i) ; % [mol m^-3 Pa^-1]
    Hs_cp_CO2(:,i) = (densf(:,i) ./ (molm_f*1e-3)) / Hv_px_Pa(4,i) ; % [mol m^-3 Pa^-1]
end

%% Calculate the molar fractions

xg_N2 = zeros(size(densf)) ;
xg_O2 = zeros(size(densf)) ;
xg_Ar = zeros(size(densf)) ;
xg_CO2 = zeros(size(densf)) ;
for i = 1 : NP
    for j = 1 : NT
        xg_N2(i,j) = (Pg_N2(i,j)*1e3) / Hv_px_Pa(1,j) ; % []
        xg_O2(i,j) = (Pg_O2(i,j)*1e3) / Hv_px_Pa(2,j) ; % []
        xg_Ar(i,j) = (Pg_Ar(i,j)*1e3) / Hv_px_Pa(3,j) ; % []
        xg_CO2(i,j) = (Pg_CO2(i,j)*1e3) / Hv_px_Pa(4,j) ; % []
    end
end

%% Calculate the molarity

Molg_N2 = zeros(size(densf)) ;
Molg_O2 = zeros(size(densf)) ;
Molg_Ar = zeros(size(densf)) ;
Molg_CO2 = zeros(size(densf)) ;
for j = 1 : NT
    for i = 1 : NP
        Molg_N2(i,j) = (Pg_N2(i,j)*1e3) .* Hs_cp_N2(i,j) ; % [mol m^-3]
        Molg_O2(i,j) = (Pg_O2(i,j)*1e3) .* Hs_cp_O2(i,j) ; % [mol m^-3]
        Molg_Ar(i,j) = (Pg_Ar(i,j)*1e3) .* Hs_cp_Ar(i,j) ; % [mol m^-3]
        Molg_CO2(i,j) = (Pg_CO2(i,j)*1e3) .* Hs_cp_CO2(i,j) ; % [mol m^-3]
    end
end

%% Calculate the mass fractions

cg_N2 = xg_N2 * (molm(1)/molm_f) * 1e6 ; % [mgrams of solute / kg of solvent]
cg_O2 = xg_O2 * (molm(2)/molm_f) * 1e6 ; % [mgrams of solute / kg of solvent]
cg_Ar = xg_Ar * (molm(3)/molm_f) * 1e6 ; % [mgrams of solute / kg of solvent]
cg_CO2 = xg_CO2 * (molm(4)/molm_f) * 1e6 ; % [mgrams of solute / kg of solvent]

% air
cg = cg_N2 + cg_O2 + cg_Ar ; % [mgrams of air / kg of water]

c_N2_mgL = cg_N2 .* densf / 1000 ; % [mgrams of solute / Litre of solvent]
c_O2_mgL = cg_O2 .* densf / 1000 ; % [mgrams of solute / Litre of solvent]
c_Ar_mgL = cg_Ar .* densf / 1000 ; % [mgrams of solute / Litre of solvent]
c_CO2_mgL = cg_CO2 .* densf / 1000 ; % [mgrams of solute / Litre of solvent]

c_air_mgL = c_N2_mgL + c_O2_mgL + c_Ar_mgL + c_CO2_mgL ; % [mgrams of air / Litre of water]

% out = [P_atm, c_N2_mgL, c_O2_mgL, c_Ar_mgL, c_CO2_mgL, c_air_mgL] ;
% out = [c_O2_mgL, c_air_mgL] ;

%% Calculate the partial molar volumes

volpm_N2 = (34.5 -0.03 .* Temp') .* 1e-6 ; % [m^3 of gas / mol of gas]
volpm_O2 = (31.7 -0.04 .* Temp') .* 1e-6 ; % [m^3 of gas / mol of gas]
volpm_Ar = (32.7 -0.06 .* Temp') .* 1e-6 ; % [m^3 of gas / mol of gas]
volpm_CO2 = 34.2 .* 1e-6 ; % [m^3 of gas / mol of gas]

%% Calculate the volumes

% molar fraction of water
xw = 1 - (xg_N2 + xg_O2 + xg_Ar + xg_CO2 ) ; % []
% molar volume of water
volm_w = (molm_f .* 1e-3) ./densf ; % [m^3 of water / mol of water]

% total molar volume of the mixture
vol_t = xw .* volm_w + xg_N2 .* volpm_N2 + ...
    xg_O2 .* volpm_O2 + xg_Ar .* volpm_Ar + xg_CO2 .* volpm_CO2 ; % [m^3 / mol]

% total molar mass of the mixture
molm_t = xw .* molm_f + xg_N2 .* molm(1) + xg_O2 .* molm(2) + ...
    xg_Ar .* molm(3) + xg_CO2 .* molm(4) ; % [gr / mol]

%% Calculate the variation in the density

delta_rho = (molm_t .* 1e-3) ./ vol_t...
    - densf ; % [kg/m^3]

%% Calculate the variation in the density due to each gas component
% [1985] H. Watanabe and K. Iizuka: The Influence of Dissolved Gases on the Density of Water

volpm_N2_B = (35.359 + 0.0075 .* Temp') .* 1e-6 ; % [m^3 of gas / mol of gas]
volpm_O2_B = (31.265 -0.0285 .* Temp') .* 1e-6 ; % [m^3 of gas / mol of gas]
volpm_Ar_B = (30.092 -0.0277 .* Temp') .* 1e-6 ; % [m^3 of gas / mol of gas]
volpm_CO2_B = (33.520 + 0.0071 .* Temp') .* 1e-6 ; % [m^3 of gas / mol of gas] look at the erratum paper, which suggests 33.520 is correct NOT 35.520

% standar volume of gas
Vs = 22.41383 .* 1e-3 ; % [m^3 mol^-1] of gas

% volume concentrations of gasses in water
ax = Vs .* (densf ./ (molm_f.*1e-3) ) ;
volc_N2 =  ax .* xg_N2 ; % [m^3 of gas / m^3 of water]
volc_O2 = ax .* xg_O2 ; % [m^3 of gas / m^3 of water]
volc_Ar = ax .* xg_Ar ; % [m^3 of gas / m^3 of water]
volc_CO2 = ax .* xg_CO2 ; % [m^3 of gas / m^3 of water]

delta_rho_N2 = (volc_N2 ./ Vs) .* ( molm(1).*1e-3 - densf .* volpm_N2_B ) ; % [kg m^-3]
delta_rho_O2 = (volc_O2 ./ Vs) .* ( molm(2).*1e-3 - densf .* volpm_O2_B ) ; % [kg m^-3]
delta_rho_Ar = (volc_Ar ./ Vs) .* ( molm(3).*1e-3 - densf .* volpm_Ar_B ) ; % [kg m^-3]
delta_rho_CO2 = (volc_CO2 ./ Vs) .* ( molm(4).*1e-3 - densf .* volpm_CO2_B ) ; % [kg m^-3]

delta_rho_B = delta_rho_N2 + delta_rho_O2 + delta_rho_Ar + delta_rho_CO2 ;  % [kg m^-3]

%% Calculate the density difference ratio percentage
rho_ratio = (delta_rho ./ densf) .* 100 ; % [%]
rho_ratio_B = (delta_rho_B ./ densf) .* 100 ; % [%]

%% OUTPUT: density changes

dens_var = [Temp, delta_rho_N2'.*1e3, delta_rho_O2'.*1e3, delta_rho_Ar'.*1e3, delta_rho_CO2'.*1e3, delta_rho_B'.*1e3, ...
    delta_rho'.*1e3, ...
    delta_rho_B'.*1e5 ./ densf' , delta_rho'.*1e5 ./ densf'] ;