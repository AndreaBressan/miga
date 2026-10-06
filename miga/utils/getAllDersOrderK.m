function allDers = getAllDersOrderK(nVars,k)
    ranges = repmat({0:k}, 1, nVars);
    [grids{1:nVars}] = ndgrid(ranges{:});
    allDers = cell2mat(cellfun(@(x) x(:).', grids, 'UniformOutput', false).');
    totalOrder = sum(allDers, 1);
    allDers = allDers(:, totalOrder <= k);
    [~, idx] = sort(sum(allDers, 1));
    allDers = allDers(:, idx);
end