function [resBasis, idx] = tRestrict(Basis,side)
arguments
    Basis (1,1) aBasis
    side  (:,1) double {mustBeInRange(side,-1,1),mustBeInteger(side)}
end
internal = find(side==0);
extremes = setdiff(1:length(side),internal);
 
if  isa(Basis,'exprBasis')
    idx=1:Basis.dimension;
   newStr = Basis.expression;

    for i = extremes
         xx = strcat('x',num2str(i));

         if side(i)==1
            new = Basis.domain(i,2);
         elseif side(i)==-1
             new = Basis.domain(i,1);
         end

       newStr = replace(newStr,xx,num2str(new));
    end

    for i = 1:length(internal)
        xx = strcat('x',num2str(internal(i)));
        new =  strcat('x',num2str(i));
        newStr = replace(newStr,xx,num2str(new));
    end


    resBasis = exprBasis(Basis.domain(internal,:),newStr).setParameters(Basis.parametersNames,Basis.parametersValue);
    

elseif  isa(Basis,'tBasis')
    idx=tSideIndices(Basis,side);

    for i=length(internal):-1:1
        deg_new = Basis.components(internal(i)).degree;
        knots_new = Basis.components(internal(i)).knots;
        periodic_new = Basis.components(internal(i)).periodic;
        domain_new = Basis.components(internal(i)).domain;
        resNew(i)=bSpline(deg_new, knots_new, periodic_new, domain_new);
    end

    if isscalar(internal)
        resBasis = resNew(1);
    else
        resBasis = tBasis(resNew);
    end
    % Sto supponendo di aver un open knot vector e che quindi il valore della prima e ultima funzione di base siano 1
   
 
    % coeff=1;
    % for idx = 1:length(extremes)
    %      if side(idx)==1
    %        coeff = coeff*Basis.components(extremes(idx)).evaluate(Basis.domain(extremes(idx),2));
    %        coeff=coeff(1,1,end);
    %      elseif side(idx)==-1
    %        coeff = coeff*Basis.components(extremes(idx)).evaluate(Basis.domain(extremes(idx),1));
    %        coeff = coeff(1,1,1);
    %      end
    % end
end
end

