x=imread('canny.jpg');
Iblur = imgaussfilt(x,2);
i=rgb2gray(x);
[r,c,d]=size(i);
i=im2double(i);
for m=2:r-1
    for n=2:c-1
        for o=1:d
            a(m,n,o)=(0)*i(m-1,n-1,o)+(-1)*i(m-1,n,o)+(0)*i(m-1,n+1,o)+(-1)*i(m,n-1,o)+(4)*i(m,n,o)+(-1)*i(m,n+1,o)+(0)*i(m+1,n-1,o)+(-1)*i(m+1,n,o)+(0)*i(m,n+1,o);
           
        end
    end
end
imshow(a);