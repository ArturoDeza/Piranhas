function im_Lab=RGB2Lab(im)
% im_Lab=RGB2Lab(im)
% converts RGB color space to CIELab color space
% im is an input RGB image
% im_Lab is the output Lab image

im=double(im); %Convert unsigned 8-bit integer to double precision.
               
r=im(:,:,1)/255; %get r,g,b value in the range of [0,1]
g=im(:,:,2)/255;
b=im(:,:,3)/255;

% the figure from graybar.m and the infromation from the website 
% http://www.cinenet.net/~spitzak/conversion/whysrgb.html, we can conclude
% that our RGB system is sRGB

% if RGB system is sRGB
indr=r>0.04045;
r(indr)=((r(indr)+0.055)/1.055).^2.4;
r(~indr)=r(~indr)/12.92;

indg=g>0.04045;
g(indg)=((g(indg)+0.055)/1.055).^2.4;
g(~indg)=g(~indg)/12.92;

indb=b>0.04045;
b(indb)=((b(indb)+0.055)/1.055).^2.4;
b(~indb)=b(~indb)/12.92;

%Observer. = 2°, Illuminant = D65
X = r * 0.4124 + g * 0.3576 + b * 0.1805;
Y = r * 0.2126 + g * 0.7152 + b * 0.0722;
Z = r * 0.0193 + g * 0.1192 + b * 0.9505;

x = X /  95.047; %x=X/Xn         %Observer = 2°, Illuminant = D65
y = Y / 100.000; %y=Y/Yn  
z = Z / 108.883; %z=Z/Zn  

indx=x>0.008856;
x(indx)=x(indx).^(1/3);
x(~indx)=7.787*x(~indx)+16/116;

indy=y>0.008856;
y(indy)=y(indy).^(1/3);
y(~indy)=7.787*y(~indy)+16/116;

indz=z>0.008856;
z(indz)=z(indz).^(1/3);
z(~indz)=7.787*z(~indz)+16/116;

CIE_L = ( 116 * y ) - 16;
CIE_a = 500 * ( x - y );
CIE_b = 200 * ( y - z );

im_Lab(:,:,1)=CIE_L;
im_Lab(:,:,2)=CIE_a;
im_Lab(:,:,3)=CIE_b;