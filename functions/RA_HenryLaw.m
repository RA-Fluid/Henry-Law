% function to find Henry's volatility law Constant
% inputs are the gas name and temperature in kelvin
% Output: in Pa
% source: 
% T. R. Rettich, R. Battino, and E. Wilhelm, J. Solution Chem. 13:335 (1984)
% T. R. Rettich, R. Battino, and E. Wilhelm, J. Chem. Thermodyn. 32:1145 (2000)
% T. R. Rettich, R. Battino, and E. Wilhelm, J. Solution Chem. 21:987 (1992)
% J. J. Carroll, J. D. Slupsky, and A. E. Mather, J. Phys. Chem. Ref. Data 20:1201 (1991)

function Hv_px = RA_HenryLaw(gnam,T)

if strcmp(gnam,'carbon dioxide') == 1

    a_0 = 6.9809;
    a_1 = 1.2817e4; % [K]
    a_2 = -3.7668e6; % [K^2]
    a_3 = 2.997e8; % [K^3]

elseif strcmp(gnam,'nitrogen') == 1

    a_0 = 14.2766192;
    a_1 = 6.3866654e3;
    a_2 = -1.1397892e6;
    a_3 = 0;

elseif strcmp(gnam,'oxygen') == 1

    a_0 = 14.989460;
    a_1 = 5.742622e3;
    a_2 = -1.070683e6;
    a_3 = 0;

elseif strcmp(gnam,'Argon') == 1

    a_0 = 15.349542;
    a_1 = 5.467601e3;
    a_2 = -1.029186e6;
    a_3 = 0;

end

Hv_px = exp( a_0 + a_1/T + a_2/(T^2) + a_3/(T^3) ) ; % [Pa]