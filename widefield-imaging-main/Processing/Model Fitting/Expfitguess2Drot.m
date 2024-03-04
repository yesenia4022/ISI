
function G = Expfitguess2Drot

%Double check Initial guesses

global RF origuess

sigguess = length(RF(:,1))/5;

Base = prctile(RF(:),10);
Amp = max(RF(:)-Base);

G = [sigguess sigguess Amp Base origuess];


%%%%%
% imW = length(RF(:,1));
% 
% [x y] = meshgrid(domx,domy);
% 
% xp = x*cos(G(5)*pi/180) + y*sin(G(5)*pi/180);
% yp = y*cos(G(5)*pi/180) - x*sin(G(5)*pi/180);
% 
% ffitD1 = exp(-(yp-G(1)).^2/(2*G(3).^2));
% ffitD2 = exp(-(xp-G(2)).^2/(2*G(4).^2));
% ffit = ffitD1.*ffitD2;
% 
% figure,imagesc(RF),colorbar
% figure,imagesc(ffit),colorbar
% err = mean((ffit(:)-RF(:)).*(ffit(:)-RF(:)))
% %%%%