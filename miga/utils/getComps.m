function comps = getComps(delta)

elements = 1:length(delta);
seed = repelem(elements, delta);
comps = unique(perms(seed), 'rows')';

end