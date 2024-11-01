% Function to create an array with 5 climate types from the
% Koppen Geiger climate classes (Beck et al. 2023)
% 1. Humid tropics (Af, Am)
% 2. Wet-dry tropics (Aw)
% 3. (Semi)arid (sub)tropics (BWh, BWk, BSh, BSk)
% 4. Humid subtropics (Cwa, Cfa)
% 5. Highland (sub)tropics (Cwb, Cfb, Dwb, Dwc, EF, ET)

function[T]=climate_classification(T)
T.Classes=T.Climate;
T.Classes(T.Classes==22)=0;
T.Classes(T.Classes==23)=0;
T.Classes(T.Classes==29)=0;
T.Classes(T.Classes==30)=0;
T.Classes(T.Classes==9)=NaN;
T.Classes(T.Classes>15)=NaN;
T.Classes(T.Classes==1)=1;
T.Classes(T.Classes==2)=1;
T.Classes(T.Classes==3)=2;
T.Classes(T.Classes==4)=3;
T.Classes(T.Classes==5)=3;
T.Classes(T.Classes==6)=3;
T.Classes(T.Classes==7)=3;
T.Classes(T.Classes==11)=4;
T.Classes(T.Classes==14)=4;
T.Classes(T.Classes==12)=5;
T.Classes(T.Classes==15)=5;
T.Classes(T.Classes==0)=5;
T.Classes=categorical(T.Classes);
T.KG=T.Classes;

clear T.Classes

end
