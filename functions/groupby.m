function groups = groupby(strs)
    
    gs = cellfun(@(x) find(strcmp(x, strs)), unique(strs), 'UniformOutput', false);
    groups = [unique(strs), gs];
    
end