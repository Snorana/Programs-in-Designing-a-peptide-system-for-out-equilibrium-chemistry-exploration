function [comp2,Z]=aacharge(monoseq,comp,forme,pH)
%fonction qui prend en argument un tableau de composition en acide aminé,
%la forme linéaire ou cyclique et le pH de la solution
[l,~] = size(comp);
comp2 = cell(l+10,2);
comp2(1:l,:)=comp(:,:);
ligne = 0;
Z=0;
r = 0;
b  = 0;
aaseq = num2cell(monoseq);
c = length(aaseq);
backbone = zeros(1,2);
%on initialise le tableau contenant l'acide aminé, peut il être chargé ou
%non le pKa de l'acide, de l'amine, du résidu et la charge du résidu après
%réaction acide/base.
SCpK={'A','o',2.34,9.69,0,0;
    'R','R+',2.17,9.04,12.48,1;
    'N','o',2.02,8.80,0,0;
    'D','D-',1.88,9.60,3.65,-1;
    'C','C-',1.96,10.28,8.18,-1;
    'E','E-',2.19,9.67,4.25,-1;
    'Q','o',2.17,9.13,0,0;
    'G','o',2.34,9.6,0,0;
    'H','H+',1.82,9.17,6,1;
    'I','o',2.36,9.6,0,0;
    'L','o',2.36,9.6,0,0;
    'K','K+',2.18,8.95,10.53,1;
    'M','o',2.28,9.21,0,0;
    'F','o',1.83,9.13,0,0;
    'P','o',1.99,10.6,0,0;
    'S','o',2.21,9.15,0,0;
    'T','o',2.09,9.1,0,0;
    'W','o',2.83,9.39,0,0;
    'Y','Y-',2.2,9.11,10.07,-1;
    'V','o',2.32,9.62,0,0};
[aa,~] = size(SCpK);
%pour chaque aide aminé du tableau de composition on l'identifie dans le
%tableau ci-dessus et on calcul sa quantité sous forme acide et base et on
%crée une ligne avec la quantité chargée en fin de tableau.
for i = 1:l
    for j=1:aa
        if strcmp(comp(i,1),SCpK(j,1))==1 && strcmp(SCpK(j,2),'o') ~=1 
            a = 1/(1+10^(cell2mat(SCpK(j,6))*(pH-cell2mat(SCpK(j,5)))));
            Z= Z+ cell2mat(comp(i,2))*a*cell2mat(SCpK(j,6));
            ligne = ligne + 1;
            comp2(ligne+l,1) = SCpK(j,2);
            comp2{ligne+l,2} = cell2mat(comp(i,2))*a;
            comp2{i,2} = cell2mat(comp(i,2))*(1-a);
            %on identifie le pKa de l'amine et de l'acide carboxilique
            %terminaux
            if r ==0 || b ==0
                if strcmp(aaseq(1),SCpK(j,1))==1
                    backbone(2) = cell2mat(SCpK(j,5));
                    r=1;
                end
                if strcmp(aaseq(c),SCpK(j,1))==1
                    backbone(1) = cell2mat(SCpK(j,3));
                    b=1;
                end
            end
        end
    end
end
%reste à calculer sous quel forme se trouve l'amine et l'acide carboxilique
%terminaux si le peptide est sous forme linéaire
if forme == 1
            if pH<(backbone(1)+backbone(2))/2
                a = 1/(1+10^((pH-backbone(1))));
                ligne = ligne+1;
                comp2{ligne+l,1} = 'AABB';
                comp2{ligne+l,2} = 1-a;
                ligne = ligne+1;
                comp2{ligne+l,1} = 'AABB+';
                comp2{ligne+l,2} = a;
            else
                a = 1/(1+10^(-(pH-backbone(2))));
                ligne = ligne+1;
                comp2{ligne+l,1} = 'AABB';
                comp2{ligne+l,2} = 1-a;
                ligne = ligne+1;
                comp2{ligne+l,1} = 'AABB-';
                comp2{ligne+l,2} = a;
            end
end
comp2 = comp2(1:l+ligne,:);
            
            
    
