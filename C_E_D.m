%%canny edge detection
%%step - 1 :-> smoothing using gaussian filter
clear all;
close all;
clc;
A = imread('remote.jpeg');

B = r2g(A);

imshow(B); title('Original GrayScale Image'); figure();
B = double(B);

[r,c] = size(B);

r1 = 2*r;
c1 = 2*c;

p = zeros(r1,c1);
q = zeros(r1,c1);

for i = 1:r
    for j = 1:c
        p(i,j) = B(i,j);
    end
end

for i = 1:r1
    for j = 1:c1
        q(i,j) = (-1)^(i+j)*p(i,j);
    end
end

f = fft2(q);
thresh = 200;
t= glp(f,thresh);
u = ifft2(f);

%%remove the padding

for i = 1:r
    for j = 1:c
        s(i,j) = u(i,j);
    end
end

s = uint8(abs(s));
imshow(s); title('GAUSSIAN FILTERED IMAGE');
figure();
 
%%step-2 :-> sobel edge detector
s = double(s);
for i = 1:r-2
    for j = 1:c-2
        %sobel mask for x direction:
        Gx=((2*s(i+2,j+1)+s(i+2,j)+s(i+2,j+2))-(2*s(i,j+1)+s(i,j)+s(i,j+2)));
        %Sobel mask for y-direction:
        Gy=((2*s(i+1,j+2)+s(i,j+2)+s(i+2,j+2))-(2*s(i+1,j)+s(i,j)+s(i+2,j)));
s(i,j)=sqrt(Gx.^2+Gy.^2);
s_mag(i,j) = sqrt(Gx.^2+Gy.^2);
n_direction(i,j) = atan2(Gy,Gx);
 n_direction(i,j) = n_direction(i,j)*(180/pi);
    end
end

imshow(s);
figure();
%threshold =  100;
%Define a threshold value
Thresh1=120;
s=max(s,Thresh1);
s(s==round(Thresh1))=0;
s=uint8(s);
imshow(~s); figure(); title('Edge detected Image');
%n_direction = zeros(r-2,c-2);
%%discretisation of direction
%for i = 1:r-2
 %   for j = 1:c-2
  %     n_direction = atan2(Gy,Gx)
   %    n_direction = n_direction*(180/pi);
    %end 
%end
       n_direction_dis = zeros(r-2,c-2);
for i = 1:r-2
    for j = 1:c-2
        if ((n_direction(i, j) > 0 ) && (n_direction(i, j) < 22.5) || (n_direction(i, j) > 157.5) && (n_direction(i, j) < -157.5))
            n_direction_dis(i, j) = 0;
        end
        
         if ((n_direction(i, j) > 22.5) && (n_direction(i, j) < 67.5) || (n_direction(i, j) < -112.5) && (n_direction(i, j) > -157.5))
            n_direction_dis(i, j) = 45;
        end
        
        if ((n_direction(i, j) > 67.5 && n_direction(i, j) < 112.5) || (n_direction(i, j) < -67.5 && n_direction(i, j) > 112.5))
            n_direction_dis(i, j) = 90;
        end
        
        if ((n_direction(i, j) > 112.5 && n_direction(i, j) <= 157.5) || (n_direction(i, j) < -22.5 && n_direction(i, j) > -67.5))
            n_direction_dis(i, j) = 135;
        end
    end
end
imagesc(n_direction_dis); colorbar;
title('Normal Directions');

%%step-4 :->non_maximum suppression
sup_im = zeros(r-2,c-2);
for i = 2:r-3
    for j = 2:c-3
%         if(n_direction(i,j)==0)
% if(((s_mag(i,j)>s_mag(i,j-1)||(s_mag(i,j)>s_mag(i,j+1)))
%     sup_im(i,j) = s_mag(i,j);
% else
%     sup_im(i,j) =0;
% end
% end
%         
%  if(n_direction(i,j)==45)        
%  if((s_mag(i,j)>s_mag(i+1,j-1)||(s_mag(i,j)>(s_mag(i-1,j+1))
%     sup_im(i,j) = s_mag(i,j);
% else
%     sup_im(i,j) =0;
% end
% end
%         
% if(n_direction(i,j)==90)
% if((s_mag(i,j)>s_mag(i-1,j)||(s_mag(i,j)>(s_mag(i+1,j))
%     sup_im(i,j) = s_mag(i,j);
% else
%     sup_im(i,j) =0;
% end
% end
% 
% if(n_direction(i,j)==135)
% if((s_mag(i,j)>s_mag(i-1,j-1)||(s_mag(i,j)>(s_mag(i+1,j+1))
%     sup_im(i,j) = s_mag(i,j);
% else
%     sup_im(i,j) =0;
% end
% end         
% 
%     end
% end
% figure();
% imshow(sup_im);
 if (n_direction_dis(i, j) == 0)
        	if (s_mag(i, j) > s_mag(i, j - 1) && s_mag(i, j) > s_mag(i, j + 1))
                sup_im(i, j) = s_mag(i, j);
            else
                sup_im(i, j) = 0;
            end
        end

        if (n_direction_dis(i, j) == 45)
            if (s_mag(i, j) > s_mag(i + 1, j - 1) && s_mag(i, j) > s_mag(i - 1, j + 1))
                sup_im(i, j) = s_mag(i, j);
            else
                sup_im(i, j) = 0;
            end
        end

        if (n_direction_dis(i, j) == 90)
            if (s_mag(i, j) > s_mag(i - 1, j) && s_mag(i, j) > s_mag(i + 1, j))
                sup_im(i, j) = s_mag(i, j);
            else
                sup_im(i, j) = 0;
            end
        end

        if (n_direction_dis(i, j) == 135)
            if (s_mag(i, j) > s_mag(i - 1, j - 1) && s_mag(i, j) > s_mag(i + 1, j + 1))
                sup_im(i, j) = s_mag(i, j);
            else
                sup_im(i, j) = 0;
            end
        end

    end
end

figure();
imshow(sup_im); title('Suppressed image');

%%step - 5 :-> hysterisis thresholding
L_T_F = 0.1;
H_T_F = 0.5;
ThreshL = L_T_F * max(max(sup_im));
ThreshH = H_T_F * max(max(sup_im));

thresh_im = zeros(r, c);
for i = 2  : r-2
    for j = 2 : r-2
        if (sup_im(i, j) < ThreshL)
            thresh_im(i, j) = 0;
        elseif (sup_im(i, j) > ThreshH)
            thresh_im(i, j) = 1;
        else
            if (  sup_im(i,j)>ThreshH &&                                     ...
                 ((sup_im(i + 1, j) > ThreshL)&&(sup_im(i + 1, j) < ThreshH) )|| ...
                 ((sup_im(i + 1, j) > ThreshL) &&(sup_im(i + 1, j+1) < ThreshH)) || ...
                 ((sup_im(i + 1, j) > ThreshL)&&(sup_im(i + 1, j-1) < ThreshH)) || ...
                 ((sup_im(i - 1, j-1) > ThreshL)&&(sup_im(i + 1, j) < ThreshH))  || ...
                 ((sup_im(i - 1, j+1) > ThreshL)&&(sup_im(i + 1, j) < ThreshH) ) || ...
                 ((sup_im(i - 1, j) > ThreshL)&&(sup_im(i + 1, j) < ThreshH))  || ...
                 ((sup_im(i, j + 1) > ThreshL)&&(sup_im(i + 1, j) < ThreshH))  || ...
                 ((sup_im(i, j - 1) > ThreshL)&&(sup_im(i + 1, j) < ThreshH))     ...
               )
                thresh_im(i-1, j) = 1;
                thresh_im(i+1, j) = 1;
                thresh_im(i-1, j+1) = 1;
                thresh_im(i-1, j-1) = 1;
                thresh_im(i, j-1) = 1;
                thresh_im(i, j+1) = 1;
                thresh_im(i+1, j+1) = 1;
                thresh_im(i+1, j-1) = 1;
                thresh_im(i,j) = 1;
           %s(i,j) = thresh_im(i,j);
            end
        end
    end
    s(i,j) = thresh_im(i,j);
end

figure();
imshow(thresh_im);
title('Hysteresis Thresholding');
figure();
Thresh = 60;
s=max(s,Thresh);
s(s==round(Thresh))=0;

%B=uint8(B);
s = uint8(s);
imshow(s);%% title('CANNY EDGE DETECTION')