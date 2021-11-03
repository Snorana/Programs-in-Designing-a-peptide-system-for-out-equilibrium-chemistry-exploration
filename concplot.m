function concplot(A0,B0,k0,k1,tmax,Cdata)
[l,c] = size(Cdata);
d = (A0+B0)^2-4*A0*B0*(1-k1/k0);
Ce = ((A0+B0)-sqrt(d))/(2*(1-k1/k0));
p=-(k0*(B0-Ce+A0-Ce)+k1*(Ce+Ce));
q = k1-k0;
for t = 1:tmax
    E=p*Ce/((p+q*Ce)*exp(-p*t)-q*Ce);
    C = Ce-E;
    plot(t,C,'.')
    hold on
end
for i = 1:l
    plot(Cdata(i,1),Cdata(i,2),'O')
end
hold off
end
    