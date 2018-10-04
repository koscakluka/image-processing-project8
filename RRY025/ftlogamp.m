function y=ftlogamp(x);
y = fftshift(log( 1 + abs(fft2(x)) ));
figure;
subplot(2,1,1);
imshow(x,[]);
title('Input Image');
subplot(2,1,2);
imshow(y,[]);
title('Log Amp of Centred DFT');
