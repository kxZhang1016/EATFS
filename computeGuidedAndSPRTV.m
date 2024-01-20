function [G_prime, SPRTV] = computeGuidedAndSPRTV(I, mRTV, label, k)
    r = floor(k / 2);
    dimX = size(I, 1);     % dimension of I in x
    dimY = size(I, 2);     % dimension of I in y
    count = zeros(dimX, dimY);
    B = zeros(dimX, dimY);
    G = zeros(dimX, dimY);
    SPRTV = zeros(dimX, dimY);
    RTV = zeros(dimX, dimY);
    mRTV_min = zeros(size(mRTV));
    sigma_alpha = 1*k;
    eps = 10e-9;

    Ixy = imgradient(I, 'sobel');

    parfor i = 1 : dimX
        for j = 1 : dimY 
            minX = max(i-r, 1);  maxX = min(i+r, dimX);
            minY = max(j-r, 1);  maxY = min(j+r, dimY);

            I_patch = I(minX:maxX, minY:maxY);
            Ixy_patch = Ixy(minX:maxX, minY:maxY);
            mRTV_patch = mRTV(minX:maxX, minY:maxY);

            % Compute the average itensity
            B(i, j) = sum(I_patch(:)) / numel(I_patch);

            mRTV_min(i, j) = min(mRTV_patch(:));
            [row, col] = find(mRTV_patch == mRTV_min(i, j), 1);

            % patch shift after select
            minX_select = max(minX+row-1-r, 1);  maxX_select = min(minX+row-1+r, dimX);
            minY_select = max(minY+col-1-r, 1);  maxY_select = min(minY+col-1+r, dimY);  

            lengthX = maxX_select - minX_select + 1;
            lengthY = maxY_select - minY_select + 1;
            isInSP = zeros(lengthX, lengthY);
            m = 0; 
            for x = minX_select : maxX_select
                m = m+1;
                n = 0;
                for y = minY_select : maxY_select
                    n = n+1;
                    if label(x, y) == label(i, j)
                        isInSP(m, n) = 1;
                        count(i, j) =count(i, j) + 1;
                    end
                end
                %clear n;
            end
            %clear m;

            I_window = I(minX_select:maxX_select, minY_select:maxY_select);
            I_windows = I_window .* isInSP;
            Ixy_window = Ixy(minX_select:maxX_select, minY_select:maxY_select);
            Ixy_windows = Ixy_window .* isInSP;

            G(i,j) = sum(I_windows(:)) / count(i, j);

%             SPRTV(i, j) = (max(I_windows(:)) - min(I_windows(:))) * (sum(Ixy_windows(:)) / count(i, j));
            SPRTV(i, j) = (max(I_windows(:)) - min(I_windows(:))) * (sum(Ixy_windows(:)) / (k * k));
            RTV(i, j) = (max(I_patch(:)) - min(I_patch(:))) * (sum(Ixy_patch(:)) / (k * k));
            %clear isInSP;
        end
    end
    alpha = 2 * ((1 ./ (1 + exp(-sigma_alpha * (RTV - SPRTV)))) - 0.5);  
    
    % Compute the final guidance image
    G_prime = zeros(size(G));
    
    for i = 1 : size(G, 3)
        G_prime(:, :, i) = alpha .* G(:, :, i) + ...
            (1 - alpha) .* B(:, :, i);
    end
end