stim = randn(16, 16, 10000);

[spikes, gab] = Model_Neuron_V1( stim );


%%
[row, col, frames] = size(stim);
stim_reshaped = reshape(stim, row*col, frames);
delay = 0;
alpha = 0;

sta = STAcorr( stim_reshaped, spikes', alpha, delay );

sta = reshape(sta, row, col);

%%
figure;
subplot(1,2,1)
imagesc(gab)
axis off
axis equal
title('model RF')
subplot(1,2,2)
imagesc(sta)
axis off
axis equal
title('STA')
colormap gray
