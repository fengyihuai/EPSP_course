function wpeak = waveMmax(dwc, points)
% Finding the modulus maxima and its position in 
% Wavelet Transform Cofficients

% to store modulus maxima's location
ddw=zeros(size(dwc));
% to store positive maxima
pddw=ddw;
% to store negative maxima
nddw=ddw;

% Two order difference to find the positive maxima's locations
% and value
posw=dwc.*(dwc>0);
pdw=((posw(:,1:points-1)-posw(:,2:points))<0);
pddw(:,2:points-1)=((pdw(:,1:points-2)-pdw(:,2:points-1))>0);

% Two order difference to find the negative maxima's locations
% and value
negw=dwc.*(dwc<0);
ndw=((negw(:,1:points-1)-negw(:,2:points))>0);
nddw(:,2:points-1)=((ndw(:,1:points-2)-ndw(:,2:points-1))>0);

% update the locations of maxima(positive&negative)
ddw=pddw|nddw;
ddw(:,1)=1;
ddw(:,points)=1;

% generate a vector for storing the modulus maxima
wpeak=ddw.*dwc;
wpeak(:,1)=wpeak(:,1)+1e-10;
wpeak(:,points)=wpeak(:,points)+1e-10;

return;