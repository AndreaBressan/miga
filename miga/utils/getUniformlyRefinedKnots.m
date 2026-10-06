function newKnots = getUniformlyRefinedKnots(initialKnots,numLev,r)
% GETDIADICKNOTS 
% Input:
% initialKnots = initial knot vector
% numLev = number of levels
% r = in how many equal subintervals do you want to divide each knot span
% Output: 
% newKnots = cellarray with the knots of each level


arguments
    initialKnots     (1,:) double 
    numLev           (1,1) double = 1
    r                (1,1) double = 2
end

if numLev==1
    newKnots={initialKnots};
    return;
end

% build last knot vector
inSpan=r^(numLev-1)+1;
Klast=[ initialKnots(1)];
for i=2:numel( initialKnots)
    if  initialKnots(i)> initialKnots(i-1)
        ls=linspace( initialKnots(i-1), initialKnots(i),inSpan);
        Klast=[Klast,ls(2:end)];
    else
        Klast=[Klast,initialKnots(i)];
    end
end

newKnots{numLev}=Klast;
newKnots{1}= initialKnots;
% build in between levels
for l=numLev-1:-1:2
    c= [false , Klast(1:end-1)~=Klast(2:end)];
    take= ( mod(cumsum(c),r)==0) | (~ c); 
    newKnots{l}=Klast( take);
    Klast=newKnots{l};
end
end
