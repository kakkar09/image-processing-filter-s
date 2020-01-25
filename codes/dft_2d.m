x=[4,4;5,7];
[m,n]=size(x);
y=zeros(m,n);
for k=0:m-1
    for l=0:n-1
        for p=0:m-1
            for q=0:n-1
                y(k+1,l+1)=y(k+1,l+1)+x(p+1,q+1)*exp((((-i)*2*pi*k*p)/m)+(((-i)*2*pi*l*q)/n));
            end
        end    
    end
end
disp(y)
