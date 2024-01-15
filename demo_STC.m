stim = randn(16, 16, 10000);

[spikes, gab1, gab2] = Model_Neuron_V1( stim, 'complex' );


%%
[row, col, frames] = size(stim);
stim_reshaped = reshape(stim, row*col, frames);
delay = 0;
alpha = 0;

stc = STCcorr( stim_reshaped, spikes', alpha, delay );

stc1 = reshape(stc.eigenvector(:,end), row, col);
stc2 = reshape(stc.eigenvector(:,end-1), row, col);

%%
figure;
subplot(2,2,1)
imagesc(gab1)
axis off
axis equal
title('model RF')

subplot(2,2,3)
imagesc(gab2)
axis off
axis equal

subplot(2,2,2)
imagesc(stc1)
axis off
axis equal
title('STC')

subplot(2,2,4)
imagesc(stc2)
axis off
axis equal

colormap gray
