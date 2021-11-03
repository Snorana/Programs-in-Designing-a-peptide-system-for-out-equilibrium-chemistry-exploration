function comp=aacomp(monoseq,poldeg,form)
%Fonction qui prend une séquence et son degré de polymérisation, si le
%polymère est cyclique ou linéaire et renvoie un tableau avec la quantité
%de chaque acide aminé.
aaseq = num2cell(monoseq);
%un crée un tableau qui recueilllera chaque acide aminé présent et son
%nombre
c = length(aaseq);
comp = cell(c,2);
ligne = 0;
for i=1:c
    skip = 0;
    nb = 0;
    %on test si l'acide aminé n'est pas déjà présent et n'a pas déjà été
    %compté dans le tableau
    for j=1:ligne
        if strcmp(comp(j,1),aaseq(i))==1
            skip = 1;
        end
    end
    %s'il ne l'est pas
    if skip == 0
        %on compte combien de fois il apparaît dans la séquence
        for j=i:c
            if strcmp(aaseq(j),aaseq(i))==1
                nb =nb+1;
            end
        end
        %on actualise combien notre tableau d'arrivée a de lignes remplies
        %vérifier
        ligne = ligne +1;
        %on remplie la ligne avec le nombre de ce résidu en fonction du
        %degré de polymérisation
        comp(ligne,1) = aaseq(i);
        comp{ligne,2} = nb*poldeg;
        %cas particulier des formes linéaires et cycliques avec l'acide
        %aminé du début si le degré de polymérisation est superieur à 1.
        if form == 0 && i==1
            comp{ligne,2} = (nb-1)*poldeg;
        end
        if form == 1 && poldeg == 1 && i==1
            comp{ligne,2} = nb;
        end
        if form == 1 && poldeg > 1 && i==1
            comp{ligne,2} = (nb-1)*poldeg+1;
        end
    end
end
comp = comp(1:ligne,:);

                
                
    
