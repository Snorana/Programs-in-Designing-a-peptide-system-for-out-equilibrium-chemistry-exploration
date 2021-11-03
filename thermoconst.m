function [valth,Z] = thermoconst(const,monoseq,poldeg,form,pH)
%fonction mère qui prend en argument la constante à calculer, la séquence
%peptidique son dégré de polymérisation sa forme et son pH et la
%température (à venir)
comp = aacomp(monoseq,poldeg,form);
%on fait la composition de la protéine neutre
[comp,Z] = aacharge(monoseq,comp,form,pH);
%on ajuste la composition avec les charges correctes vis-à-vis du pH
[l,~] = size(comp);
col=0;
valth = 0;
nbaa = 0;
%on identifie quel constante on souhaite calculer
if strcmp(const,'G')==1
    col = 2;
end
if strcmp(const,'H')==1
    col = 3;
end
if strcmp(const,'S')==1
    col = 4;
end
%on inititalise le tableau des valeurs thermodynamiques
T = {'A',-3.94,-13.29,16.85;
    'R',31.62,-20.73,54.75;
    'R+',14.94,-35.67,60.59;
    'N',-40.61,-63.49,35.90;
    'D',-87.64,-107.13,33.24;
    'D-',-82.29,-105.33,21.36;
    'C',4.54,-4.86,24.34;
    'C-',15.95,3.77,15.03;
    'Q',-41.4,-72.65,40.43;
    'E',-88.18,-115.61,39.22;
    'E-',-82.34,-114.94,21.88;
    'G',-6.07,-5.57,17.31;
    'H',36.46,11.01,44.02;
    'H+',28.17,3.96,48.15;
    'I',2.89,-32.39,27.72;
    'L',0.67,-34.39,28.43;
    'K',5.66,-37.65,39.3;
    'K+',-8.01,-51.16,39.84;
    'M',-35.25,-59.31,40.38;
    'F',35.44,10.31,34.62;
    'P',11.31,-4.89,27.85;
    'S',-39.05,-53.21,25.26;
    'T',58.05,21.62,37.98;
    'Y',-6.92,-38.53,37.43;
    'Y-',5.5,-30.21,23.68;
    'V',-0.46,-27.21,23.70;
    'AABB',-84.87,-119.21,21.98;
    'AABB+',-88,-119.58,31.23;
    'AABB-',-71.97,-108.69,13.97};
PBB =[-21.44,-45.22,1.62];
[aa,~]=size(T);
%on calcule la contribution pour chaque élément présent dans le tableau de
%composition de la protéine
for i = 1:l
    for j=1:aa
        if strcmp(T(j,1),comp(i,1))==1
            valth = valth + cell2mat(T(j,col))*cell2mat(comp(i,2));
        end
    end
    nbaa = nbaa + cell2mat(comp(i,2));
end
if form == 1
    nbaa = round(nbaa - 1);
    valth = valth + (nbaa-1)*PBB(col-1);
end
if form == 0
    nbaa = round(nbaa);
    valth = valth + nbaa*PBB(col-1);
end
%j'affiche pour contrôler le tableau de composition de la protéine
comp
    




