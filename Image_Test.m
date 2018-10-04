%% Initialization
clc; clear; close all;
warning off;
DEBUG = 1;

%% Functions
crop = @(image, rowsN, colsN)(image(1:rowsN, 1:colsN));

%% Variables Initialization
original = imread('cameraman.tif');

motionAmmount   = 13;                       % Ammount of linear motion in Pixels
motionAngle     = 0;                        % Counterclockwise angle for motion

[originalRows, originalColumns] = size(original);   % Size of original Image

%% Fourier of Original Image
P = 2 * originalRows;
Q = 2 * originalColumns;

originalPadded = padarray(original, [P - originalRows, Q - originalColumns], 'post');

figure;
subplot(2, 2, 1);
imshow(original, []);
title('Original Image');

subplot(2, 2, 2);
imshow(originalPadded, []);
title('Original Image Padded');

originalFourier = fftshift(fft2(originalPadded));

subplot(2, 1, 2);
imshow(log(0.01 + originalFourier), []);
title('Log of Fourier of Padded Image');


%% Motion Filter
motionFilter = fspecial('motion', motionAmmount, motionAngle);
[motionRows, motionColumns] = size(motionFilter);

motionFilterPadded = padarray(motionFilter, [P - motionRows, Q - motionColumns], 'post');

figure;
subplot(1, 2, 1);
imshow(motionFilterPadded, []);
title('Padded Motion Filter');

motionFourier = fftshift(fft2(motionFilterPadded));
subplot(1, 2, 2);
imshow(motionFourier, []);
title('Fourier of Padded Motion Filter');

%% Apply Motion Filter and make Gaussian Noise

gaussianMean            = 0;        % Gaussian filter mean
gaussianVariance        = 0;     % Gaussian filter variance

motionImage = imfilter(original, motionFilter, 'conv');
motionImageFourier = fftshift(fft2(motionImage));

figure;
subplot(2, 2, 1);
imshow(log(0.01 + motionImageFourier), []);
title('Log of Fourier of Motion Image');

subplot(2, 2, 2);
imshow(motionImage, []);
title('IFT of Motion Image Cropped');

noise = double(imnoise(motionImage, 'gaussian', gaussianMean, gaussianVariance)) - double(motionImage);

subplot(2, 2, 3);
imshow(noise, []);
title('Gaussian Additive Noise');

[noiseRows, noiseColumns] = size(noise);
noiseFourier = fftshift(fft2(padarray(noise, [P - noiseRows, Q - noiseColumns], 'post')));
subplot(2, 2, 4);
imshow(log(0.01 + abs(noiseFourier)), []);
title('Fourier of Padded Noise');

%% Get Final Image
degradedImage = double(motionImage) + noise;

figure;
subplot(1, 2, 1);
imshow(degradedImage, []);
title('Image With Motion and Noise');

degradedPadded = padarray(degradedImage, [P - originalRows, Q - originalColumns], 'post');
degradedFourier = fftshift(fft2(degradedPadded));

subplot(1, 2, 2);
imshow(log(0.001 + abs(degradedFourier)), []);

%% DEBUGGING
if(DEBUG == 0)
    close all;
end
    
%% Pseudo-Inverse Filter

coefPercent = 0.7;

H = motionFourier;

H_prime = sort(abs(reshape(H, [1, numel(H)])));
pseudoThreshold = H_prime(floor(coefPercent .* numel(H)) + 1);

pseudoInverse = (abs(H) > pseudoThreshold) .* (1./H);

restorationPseudoFourier = pseudoInverse .* degradedFourier;
figure;
subplot(1, 3, 1);
imshow(original, []);
title('Original Image');

subplot(1, 3, 2);
imshow(degradedImage, []);
title('Degraded Image');

restoredPseudo = crop(ifft2(ifftshift(restorationPseudoFourier )), originalRows, originalColumns);

subplot(1, 3, 3);
imshow(restoredPseudo, []);
title('Restored Image - Pseudo Inverse');

%% Making of Wiener Filter

noisePowerSpectrum = abs(noiseFourier).^2;
originalPowerSpectrum = abs(originalFourier).^2;
K = 0;
wienerFilter = (1 ./ H) .* (abs(H).^2 ./ (abs(H).^2 + noisePowerSpectrum./originalPowerSpectrum));

restorationWienerFourier = degradedFourier .* wienerFilter;
figure;
subplot(1, 3, 1);
imshow(original, []);
title('Original Image');

subplot(1, 3, 2);
imshow(degradedImage, []);
title('Degraded Image');

restoredWiener = crop(ifft2(ifftshift(restorationWienerFourier)), originalRows, originalColumns);

subplot(1, 3, 3);
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