function allDers = getAllDers(der)
    ranges = arrayfun(@(x) 0:x,der,'UniformOutput',false);
    [grids{1:length(der)}] = ndgrid(ranges{end:-1:1});
    grids = grids(end:-1:1);
    allDers = cell2mat(cellfun(@(x) x(:).',grids,'UniformOutput',false).');
end