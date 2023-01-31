%% Converts numerical indices to logical indices (opposite of built-in find)
function idxBool = antifind(idx, len)
    idxBool = zeros(len, 1, 'logical');
    idxBool(idx) = true;
end