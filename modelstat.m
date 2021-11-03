function modelstat(C0,k1,k2,kc,tm,dt,sample,lmax)
vl = zeros(1,lmax);
vl(1) = C0;
cyst = 0;
elapsedt=cputime;
ycspecie=zeros(sample,lmax-1);
ylspecie=zeros(sample,lmax);
ylspecie(1,1)=1000*C0;
vc = zeros(1,lmax-1);
ycysteine = zeros(1,sample);
xtime = zeros(1,sample);
%vc(1:3) = C0/10;
e = 30;
p = 1;
for t=1:dt:tm
    t;
    vl2 = vl;
    vc2 = vc;
    cyst2=cyst;
    for i=1:lmax
        %Attaque des extrémités Nterminales des polymères eéquiprobable sur
        %toutes les cystéines "vulnérables" des autres polymères linéaires
        for j=1:lmax
            n = k1*vl(i)*vl(j);
            if j>1 && n>10^(-e)
                for k=1:j-1
                    if i+k<=lmax
                        vl2(i+k)= vl2(i+k) + n;
                        vl2(j-k)=vl2(j-k)+ n;
                        vl2(i) = vl2(i)-n;
                        vl2(j) = vl2(j) - n;
                    end
                end
            end
            if n>10^(-e) && i+j <= lmax
                if i==j
                    vl2(i+j)= vl2(i+j) +n/2;
                    cyst2 = cyst2 +n/2;
                    vl2(i) = vl2(i)-n/2;
                    vl2(j) = vl2(j) - n/2;
                end                
                vl2(i+j)= vl2(i+j) +n;
                cyst2 = cyst2 +n;
                vl2(i) = vl2(i)-n;
                vl2(j) = vl2(j) - n;
            end
        end
        %attaque des polymères linéaires sur les polymères cycliques, peut
        %importe la cystéine attaquée le résultat est le même
        for j=1:lmax-1
            n = k1*vl(i)*vc(j);
            if n>10^(-e) && i+j+1<= lmax
                vl2(i+j+1)= vl2(i+j+1)+ (i+1)*n;
                vc2(j) = vc2(j)-(i+1)*n;
                vl2(i) = vl2(i)-(i+1)*n;
            end
        end
        %si le polymère linéaire est plus grand que le stade monomère alors
        %ceui-ci peut cycliser
        if i>1
            n=kc(i-1)*vl(i);
            for j=2:i
            if n>10^(-e)
                vc2(j-1) = vc2(j-1) +n;
                vl2(i)= vl2(i)-n;
                if i-j>0
                vl2(i-j)=vl2(i-j)+n;
                else
                cyst2 = cyst2 +n;
                end
            end
            end
        end
    end
    %Réaction de dépolymérisation avec la nmetcyst produite avec les
    %polymères linéaires
    for i=2:lmax
        n = k2*vl(i)*cyst;
        if n>10^(-e)
            for j=1:(i-1)
                vl2(j) = vl2(j) + n;
                vl2(i-j) = vl2(i-j) + n;
                cyst2 = cyst2 - n;
                vl2(i) = vl2(i) - n;
            end
        end
    end
    %avec les polymères cycliques
    for i=1:lmax-1
        n = k2*vc(i)*cyst;
        if n>10^(-e)
            vl2(i+1)= vl2(i+1) + (i+1)*n;
            cyst2 = cyst2 - (i+1)*n;
            vc2(i) = vc2(i)-(i+1)*n;
        end
    end       
    vl=vl2;
    vc=vc2;
    cyst = cyst2;  
    if mod(t,floor(tm/sample)) == 0
        %pause()
        p=p+1;
        ylspecie(p,:)=1000*vl;
        ycspecie(p,:)=1000*vc;
        ycysteine(p)=1000*cyst;
        xtime(p) = floor(t/60);
        %figure(1)
        %xticks(0:1:10)
        %for i=1:lmax
        %plot(i,vl(i),'X');
        %hold on
        %if i>1
        %plot(i+0.5,vc(i-1),'O');
        %end
        %end
        %hold off
%        c=0;
%        for i=1:length(vl)
%            c =c+ (i+1)*vl(i);
%        end
%        for i=1:length(vc)
%            c = c+ (i+1)*vc(i);
%        end
%        c = c+cyst;
%        c*100000000
%        pause(0.05)
    end
end
ylspecie;
ycspecie;
cputime - elapsedt
figure(2)
    plot(xtime,ycysteine);
    hold on
    for i=1:lmax
    plot(xtime,ylspecie(:,i));
    hold on
    if i>1
    plot(xtime,ycspecie(:,i-1));
    end
    hold on
    end
hold off
end
