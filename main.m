clear;
clc;

addpath('SCAC');

% Load an example image 
I = im2double(imread('...'));

% Set patch size and number of iterations (listed in the image name)
k = 3;
iter = 20;
 
% Apply the bilateral texture filter
[J,SPRTV] = bilateralTextureFilter(I, k, iter);

% Display both the input and output images
figure;
set(gcf, 'Name', 'Superpixel Bilateral Texture Filtering Result');

subplot(1,2,1); imshow(I);
title('Input Image');

subplot(1,2,2); imshow(J);
title('Result of Superpixel Bilateral Texture Filtering');

% figure; imshow(SPRTV);
