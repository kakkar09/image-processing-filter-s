x=[10,-1;-2,1];
[m,n]=size(x);
y=zeros(m,n);
for k=0:m-1
    for l=0:n-1
        for p=0:m-1
            for q=0:n-1
                y(k+1,l+1)=y(k+1,l+1)+x(p+1,q+1)*cos((p+0.5)*k*(pi/m))*cos((q+0.5)*l*(pi/n));
            end
        end    
    end
end
for p=0:m-1
    for q=0:n-1
        if p==0 &&q==0
            y(p+1,q+1)=y(p+1,q+1)*((1/m)^0.5)*((1/n)^0.5)
        elseif p~=0 && q==0
             y(p+1,q+1)=y(p+1,q+1)*((2/m)^0.5)*((1/n)^0.5)
        elseif p==0 &&q~=0
            y(p+1,q+1)=y(p+1,q+1)*((1/m)^0.5)*((2/n)^0.5)
        else
            y(p+1,q+1)=y(p+1,q+1)*((2/m)^0.5)*((2/n)^0.5)
        end
    end
end    
disp(y)
