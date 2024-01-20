% Compute the mRTV 
% I: input single-channel image
% k: odd-valued patch size

function mRTV = computeMRTV(I, k)
    % Parameters
    dimX = size(I, 1);     % dimension of I in x
    dimY = size(I, 2);     % dimension of I in y
    half_k = floor(k / 2); % half of patch size
    eps = 10e-9;           % epsilon used in Eq. (4)
    
    % Compute the gradient magnitude
    Ixy = imgradient(I, 'sobel');
    
    % Initialize the matrices
    B = zeros(dimX, dimY);
    mRTV = zeros(dimX, dimY);
    
    parfor i = 1 : dimX
        for j = 1 : dimY
            
            % half_k = floor(scale(i,j) / 2);
            % Extract the local patch
            minX = max(i-half_k, 1);
            maxX = min(i+half_k, dimX);
            minY = max(j-half_k, 1);
            maxY = min(j+half_k, dimY);
            I_patch = I(minX:maxX, minY:maxY);
            Ixy_patch = Ixy(minX:maxX, minY:maxY);
            
            % Compute the mRTV based on Eq. (4)
            mRTV(i, j) = (max(I_patch(:)) - min(I_patch(:))) * ...
                max(Ixy_patch(:)) / (sum(Ixy_patch(:)) + eps);
        end
    end
end