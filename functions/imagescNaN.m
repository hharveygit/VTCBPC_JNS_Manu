% Wrapper for imagesc to plot NaN values transparent (background color)
function varargout = imagescNaN(x, y, C, varargin)
    
    bgCol = 'k';
    clim = [];
    
    % parse bgcol, clims from varargin
    ii = 1;
    while ii <= length(varargin)
        if strcmpi(varargin{ii}, 'backgroundcolor')
            bgCol = varargin{ii + 1};
            varargin(ii:ii+1) = [];
        elseif strcmpi(varargin{ii}, 'clim')
            clim = varargin{ii + 1};
            varargin(ii:ii+1) = [];
        else
            ii = ii + 2;
        end
    end

    imAlpha = ones(size(C)); % alpha mask for C
    imAlpha(isnan(C)) = 0; % set NaNs to transparent
    
    if isempty(x) || isempty(y)  
        varargout{:} = imagesc(C, 'AlphaData', imAlpha, varargin{:});
        warning('Not passing x or y inputs to imagesc()');
    else
        varargout{:} = imagesc(x, y, C, 'AlphaData', imAlpha, varargin{:});
    end
    
    if ~isempty(clim), caxis(clim); end
    set(gca, 'Color', bgCol);
    
end