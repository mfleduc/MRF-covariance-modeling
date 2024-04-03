% clear variables; close all

% addpath(genpath('NeedMat'));
% addpath(genpath('spherepts'));
% addpath(genpath('Spherical-Harmonic-Transform'));
% 
% 
% j = 0:4;
% 
% 
% theta = (0:180)*pi/180;phi = (0:359.9)*pi/180 ;
% [T,P] = ndgrid(theta, phi);
% A = get_A(2, j(1),j(end),T(:),P(:),100);
gc = zeros(size(A));
lastNdx=0;
c = 2*(0.6).^j;
dists = [];
for ii = j

    nside = get_Nside(2, ii) ;
    tp = pix2ang(nside);
    tp = cell2mat(tp);
    for kk = 1:length(tp)
        dists(:,kk) = distance([pi/2-tp(1,kk),tp(2,kk)],[pi/2-T(:),P(:)],[1 0],'radians');
    end
    gc(:,lastNdx+1:lastNdx+size(tp,2)) = GaspariCohn(c(ii+1) , dists );
    lastNdx = lastNdx+length(tp);
end

Agc = A.*gc;

save('gcneedlets_1deg_2_0_6.mat', 'Agc')