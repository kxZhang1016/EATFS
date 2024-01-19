function [J, SPRTV, mRTV_to_test, alpha] = bilateralTextureFilter(I, k, iter)
    % Verify the input image I is valid
    if ~isfloat(I) || ~sum([1, 3] == size(I, 3)) ...
            || min(I(:)) < 0 || max(I(:)) > 1
        error(['Input image I must be a double precision matrix of ', ...
               'size MxNx1 or MxNx3 on the closed interval [0,1].']);
    end
    
    % Verify patch size k
    if isempty(k) || numel(k) ~= 1 || k < 1 || mod(k, 2) ~= 1
        error('Patch size k must be at least 3 and an odd-valued.');
    end
    k = round(k);

    % Verify the number of ieration iter
    if isempty(iter) || numel(iter) ~= 1 || iter < 1
        error('There must be at least one iteration');
    end
    
    iter = round(iter);
    
    % Image parameters
    dimX = size(I, 1); % dimension of I in x
    dimY = size(I, 2); % dimension of I in y
    c = size(I, 3);    % number of color channels of I
    
    % Joint bilateral filtering parameters
    s = 2 * k - 1;             % spatial kernel size
    half_s = floor(s / 2);     % half of spatial kernel size
    sigma_s = k - 1;           % spatial sigma
    sigma_r = 0.025 * sqrt(c); % range sigma
    
    % Pre-compute the Gaussian spatial kernel
    f = fspecial('gaussian', s, sigma_s);
    
    % Initialize matrices
    B = zeros(size(I));     % blurred image
    G = zeros(size(I));
    
    SPRTVs = zeros(size(I));
    mRTVs = zeros(size(I));
    J = zeros(size(I));     % output filtered image
    
    % Superpixel Setting
%     spn=600;
    spn = floor(0.002*dimX*dimY);
    S2=1;
    ItrSet=10;
    lambda=5;
    img=uint8(I);
    [label, ~] = mex_SCAC(img, spn, S2, ItrSet, lambda);
    
    % Apply the filter 'iter' times
    for m = 1 : iter
        
        fprintf('Image Processing in Iter: %d \n',m);
        % Compute the mRTV
        for i = 1 : c
            mRTVs(:, :, i) = computeMRTV(I(:, :, i), k);
        end  
        mRTV = sum(mRTVs, 3) / c;
        
        for i = 1 : c
            [G_prime(:, :, i), SPRTVs(:, :, i)] = computeGuidedAndSPRTV(I(:, :, i), mRTV, label, k);
        end
        SPRTV = sum(SPRTVs, 3) / c;
        
        % Compute the joint bilateral filtering        
        parfor i = 1 : dimX
            for j = 1 : dimY

                minX = max(i-half_s, 1);
                maxX = min(i+half_s, dimX);
                minY = max(j-half_s, 1);
                maxY = min(j+half_s, dimY);
                
                for k = 1 : c
                    
                    % Extract the local patch
                    G_prime_patch = G_prime(minX:maxX, minY:maxY, k);
                    
                    % Compute Gaussian range kernel
                    g = exp(-(G_prime_patch - G_prime(i, j, k)).^2 ...
                            /(2 * sigma_r^2));
                    fg = g .* f((minX:maxX)-i+half_s+1, ...
                        (minY:maxY)-j+half_s+1);
                    
                    % Compute the bilateral texture filter response
                    I_patch = I(minX:maxX, minY:maxY, k);
                    J(i, j, k) = sum(I_patch(:) .* fg(:)) / sum(fg(:));
                    
                end
            end
        end
        
        % Update I for the new iteration
        I = J;
    end
end