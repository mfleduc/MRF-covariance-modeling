function y = GaspariCohn(c, d)
%GASPARICOHN Evaluation of the Gaspari-Cohn correlation function with
%cutoff radius c
% https://rmets.onlinelibrary.wiley.com/doi/10.1002/qj.49712555417
d = abs(d/c);
y = zeros(size(d));
mask1 = d<1;
mask2 = (d<=2 &d>= 1); 


y(mask1) = -0.25*(d(mask1)).^5 +0.5*(d(mask1)).^4 ...
    +5/8*(d(mask1)).^3 -5/3*(d(mask1)).^2 +1 ;
y(mask2) = 1/12*d(mask2).^5 -0.5*d(mask2).^4+5/8*d(mask2).^3 ... 
    +5/3*d(mask2).^2 ...
    -5*d(mask2)+4 -2/3*(1./d(mask2));

end

