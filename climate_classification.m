% Function to create an array with 5 climate types from the
% Koppen Geiger climate classes (Beck et al. 2023)
% 1. Humid tropics (Af, Am)
% 2. Wet-dry tropics (Aw)
% 3. (Semi)arid (sub)tropics (BWh, BWk, BSh, BSk)
% 4. Humid subtropics (Cwa, Cfa)
% 5. Highland (sub)tropics (Cwb, Cfb, Dwb, Dwc, EF, ET)

function[T]=climate_classification(T)

c=T.Climate;
c(ismember(c,[22 23 29 30]))=0;
c(ismember(c,[9,16:max(c)]))=NaN;

T.KG=NaN(size(c));
T.KG(ismember(c,[1 2]))=1;          % Humid tropics
T.KG(c==3)=2;                       % Wet-dry tropics
T.KG(ismember(c,4:7))=3;            % Arid (sub)tropics
T.KG(ismember(c,[11 14]))=4;        % Humid subtropics
T.KG(ismember(c,[0 12 13 15]))=5;   % Highland (sub)tropics

T.KG=categorical(T.KG);

end
