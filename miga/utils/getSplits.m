function splits = getSplits(alpha,n)
m = size(alpha,1);

if n == 1
    splits = reshape(alpha,m,1,1);
else
    betas = getAllDers(alpha);
    splits = zeros(m,n,0);

    for b = 1:size(betas,2)
        beta = betas(:,b);
        rest = alpha - beta;
        if any(rest < 0)
            continue;
        end
        subSplits = getSplits(rest,n-1);
        Ns = size(subSplits,3);
        if Ns == 0
            continue;
        end
        newBlock = zeros(m, n, Ns);
        newBlock(:,1,:) = repmat(beta,1,1,Ns);
        newBlock(:,2:end,:) = subSplits;
        splits(:,:,end+1:end+Ns) = newBlock;
    end
end
end
