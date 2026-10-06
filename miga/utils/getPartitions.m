function partitions = getPartitions(alpha)
    % Main entry point
    allBlocks = getAllDers(alpha); 
    % Filter out the zero block if it exists in allBlocks
    allBlocks(:, all(allBlocks == 0, 1)) = []; 
    
    partitions = recursePartitions(alpha, allBlocks, 1);
end

function partitions = recursePartitions(remaining, allBlocks, startIndex)
    if all(remaining == 0)
        partitions = {[]}; 
        return
    end

    partitions = {};
    % Only iterate from startIndex to avoid permutations (non-ordered)
    for i = startIndex:size(allBlocks, 2)
        block = allBlocks(:, i);
        
        if all(block <= remaining)
            % Recurse with the updated remaining alpha and the current index i
            subParts = recursePartitions(remaining - block, allBlocks, i);
            
            for j = 1:numel(subParts)
                partitions{end+1} = [block, subParts{j}]; %#ok<AGROW>
            end
        end
    end
end