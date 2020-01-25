i=imread('canny.jpg');

grayImage=rgb2gray(i);

gaussian1 = fspecial('Gaussian', 21, 15);

gaussian2 = fspecial('Gaussian', 21, 20);

dog = gaussian1 - gaussian2;

dogFilterImage = conv2(double(grayImage), dog, 'same');

figure;

imshow(dogFilterImage);