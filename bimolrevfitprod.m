function bimolrevfitprod(data,k0,A0,B0,r1,step,r2)
[l,~]=size(data);
minvar = -1;
min = 0;
fitk1 = 0;
for i=r1:step:r2
    k1=k0*(A0-i)*(B0-i)/i^2;
    cumvar = 0;
    for j=1:l
        p=-(k0*(B0-i+A0-i)+k1*(i+i));
        q = k1-k0;
        E=p*i/((p+q*i)*exp(-p*data(j,1))-q*i);
        Edata = i-data(j,2);
        cumvar = cumvar + (E-(i-data(j,2)))^2/E;
        %figure(1)
        %plot(data(j,1),i-E,'X')
        %hold on
        %plot(data(j,1),data(j,2),'O')
    end
    if minvar == -1
        minvar = cumvar;
    end
    if minvar>cumvar
        min = i;
        minvar = cumvar;
        fitk1 = k1;
    end
    %hold off
    %figure(2)
    plot(i,cumvar,'X');
    hold on
    %pause
end
disp(min)
disp(fitk1)
end
        
        
    
