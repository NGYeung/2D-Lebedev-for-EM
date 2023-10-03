%This script is for testing the angle
% we need the angle \in [0,\pi/2)
% V,V0,R,R0 are the orthogonal bases and sqrt(M)
function Ratio = IMRatio(V,V0,R,R0,x,y,sizes,thetaset,medium_id,fname)
%I(x,s)/I_bg(x,s)

Ratio = zeros(4,length(thetaset));
count = 1;
for theta = thetaset

a1(1) = cos(theta);
a1(2) = sin(theta);
a2(1) = -sin(theta);
a2(2) = cos(theta);

maskx = zeros(sizes.ex1(1),sizes.ex1(2));
maskx(x,y)=1;
indexx = find(maskx(:));

masky = zeros(sizes.ey2(1),sizes.ey2(2));
masky(x,y)=1;
indexy = find(masky(:))+sizes.ex1p;


Ix1 = sum(((a1(1).*V0(indexx,:)+a1(2).*V0(indexy,:))*R).^2,2);
    
Ixr1 = sum(((a1(1).*V(indexx,:)+a1(2).*V(indexy,:))*R).^2,2);
    
Ix01 = sum(((a1(1).*V0(indexx,:)+a1(2).*V0(indexy,:))*R0).^2,2);

Iy2 = sum(((a2(1).*V0(indexx,:)+a2(2).*V0(indexy,:))*R).^2,2);

Iyr2 = sum(((a2(1).*V(indexx,:)+a2(2).*V(indexy,:))*R).^2,2);

Iy02 = sum(((a2(1).*V0(indexx,:)+a2(2).*V0(indexy,:))*R0).^2,2);


R1ideal = max(Ixr1./Ix01)/min(Ixr1./Ix01);
R1rom = max(Ix1./Ix01)/min(Ix1./Ix01);
R2ideal = max(Iyr2./Iy02)/min(Iyr2./Iy02);
R2rom = max(Iy2./Iy02)/min(Iy2./Iy02);

Ratio(1,count) = R1ideal;
Ratio(2,count) = R1rom;
Ratio(3,count) = R2ideal;
Ratio(4,count) = R2rom;

count = count+1;
end

filename = strcat('results/', fname, '/RatiovAngle_', '_', medium_id,  '_', num2str(theta), '.mat');
save(filename,'Ratio')

end