function C = fdbCoeff(alpha,gamma)
remaining = alpha;
C = 1;
for k = 1:size(gamma,2)
    g = gamma(:,k);
    C = C * multiNchooseK(remaining,g);
    remaining = remaining - g;
end
[~,~,idx] = unique(gamma','rows');
counts = accumarray(idx,1);
C = C / prod(factorial(counts));
end



