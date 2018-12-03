%%RGB TO GRAY CONVERSION
GIm  = imread('CAR.jpg');
red = GIm(:,:,1);
green = GIm(:,:,2);
blue = GIm(:,:,3);
GImg = 0.21*red + 0.72*green + 0.07*blue;
imshow(GImg);

