% Group Number - p8g2

%% Initialization
clc; clear; close all;
warning off;

%% Functions
crop = @(image, rowsN, colsN)(image(1:rowsN, 1:colsN));

%% Variables Initialization
original = mat2gray(imread('cameraman.tif'));

[originalRows, originalColumns] = size(original);   % Size of original Image

%% Fourier of Original Image
P = 2 * originalRows;
Q = 2 * originalColumns;

originalPadded = padarray(original, [P - originalRows, Q - originalColumns], 'post');

figure;
subplot(1, 3, 1);
imshow(original, []);
title('Original Image');

subplot(1, 3, 2);
imshow(originalPadded, []);
title('Original Image Padded');

originalFourier = fftshift(fft2(originalPadded));

subplot(1, 3, 3);
imshow(log(0.001 + real(originalFourier)), []);
title('Log of Fourier of Original Image');

%% Motion Filter

motionAmmount   = 13;                       % Ammount of linear motion in Pixels
motionAngle     = 0;                        % Counterclockwise angle for motion

PSF = fspecial('motion', motionAmmount, motionAngle);
[motionRows, motionColumns] = size(PSF);

motionFilterPadded = padarray(PSF, [P - motionRows, Q - motionColumns], 'post');

figure;
subplot(1, 2, 1);
imshow(motionFilterPadded, []);
title('Padded Motion Filter');

H = fftshift(psf2otf(PSF, [P Q]));

subplot(1, 2, 2);
imshow(log(0.001 + H), []);
title('Fourier of Padded Motion Filter');

%% Apply Motion Filter and make Gaussian Noise

gaussianMean            = 0;     % Gaussian filter mean
gaussianDeviation       = 0;     % Gaussian filter variance

motionImage = imfilter(original, PSF, 'circular', 'conv');

% figure;
% % subplot(2, 2, 1);
% imshow(motionImage, []);
% title('Motion Image applied in Space');

% figure;
% % subplot(2, 2, 2);
% imshow(log(0.01 + motionImageFourier), []);
% title('Log of Fourier of Motion Image');

noise = randn(size(motionImage)) * gaussianDeviation + gaussianMean; % From standard distribution to chosen distribution.
[noiseRows, noiseColumns] = size(noise);
noiseFourier = fftshift(fft2(padarray(noise, [P - noiseRows, Q - noiseColumns], 'post')));

figure;
% subplot(2, 2, 3);
imshow(noise, []);
% title('Gaussian Additive Noise');

figure;
% subplot(2, 2, 4);
imshow(log(0.01 + abs(noiseFourier)), []);
title('Fourier of Padded Noise');

% %% Get Final Image
degradedImage = imadd(motionImage, noise);
[degradedRows, degradedColumns] = size(degradedImage);

figure;
subplot(1, 2, 1);
imshow(degradedImage, []);
title('Image With Motion and Noise');

degradedPadded = padarray(degradedImage, [P - degradedRows, Q - degradedColumns], 'post');
degradedFourier = fftshift(fft2(degradedPadded));

subplot(1, 2, 2);
imshow(log(0.001 + abs(degradedFourier)), []);
title('Fourier of Image With Motion and Noise');

%% DEBUGGING
if(1)
    close all;
end

%% Pseudo-Inverse Filter

coefPercent = 0.7;

H_prime = sort(abs(reshape(H, [1, numel(H)])));
pseudoThreshold = H_prime(floor(coefPercent .* numel(H)) + 1);

pseudoInverse = (abs(H) > pseudoThreshold) .* (1 ./ H);

restorationPseudoFourier = pseudoInverse .* degradedFourier;

figure;
subplot(2, 2, 1);
imshow(original, []);
title('Original Image');

subplot(2, 2, 2);
imshow(degradedImage, []);
title('Degraded Image');

restoredPseudo = crop(ifft2(ifftshift(restorationPseudoFourier )), degradedRows, degradedColumns);

subplot(2, 2, 3);
imshow(log(0.001 + restorationPseudoFourier), []);
title('Restored Fourier');

subplot(2, 2, 4);
imshow(restoredPseudo, []);
title('Restored Image - Pseudo Inverse');

%% Making of Wiener Filter

CHEATING = 0;

if(CHEATING > 0)
    noisePowerSpectrum = abs(noiseFourier).^2;
    originalPowerSpectrum = abs(originalFourier).^2;
    K = noisePowerSpectrum./originalPowerSpectrum;
else
    K = 0;
%     temp = fspecial('average', 5);
%     smoothedG = imfilter(degradedFourier, temp, 'replicate', 'same');
%     
%     [X, Y] = meshgrid(1:1:P);
%     
%     mesh(X, Y, real(smoothedG));
end
    
wienerFilter = (1 ./ H) .* (1 ./ (1 + K ./ (abs(H).^2) ));

restorationWienerFourier = degradedFourier .* wienerFilter;
figure;
subplot(2, 2, 1);
imshow(original, []);
title('Original Image');

subplot(2, 2, 2);
imshow(degradedImage, []);
title('Degraded Image');

restoredWiener = crop(ifft2(ifftshift(restorationWienerFourier)), degradedRows, degradedColumns);

subplot(2, 2, 3);
imshow(log(0.001 + restorationWienerFourier), []);
title('Restoration Fourier');

subplot(2, 2, 4);
imshow(restoredWiener, []);
title('Restored Image - Wiener');

%% Comparison

figure;
subplot(1, 2, 1);
imshow(restoredPseudo, []);
title('Pseudo Inverse Filter');
subplot(1, 2, 2);
imshow(restoredWiener, []);
title('Wiener Filter');

average = fspecial('average', 5);
averageG = imfilter(degradedFourier, average, 'conv');
