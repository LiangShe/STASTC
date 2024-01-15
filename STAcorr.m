function [ sta ] = STAcorr( stim, spk, alpha, delay )
% [ sta ] = STAcorr( stim, spk, alpha, delay )
% Corrected STA / Regularized STA 
% see ref 1 for explaination of the problem, see also ref 3
% the implementation from ref 2
% 
% INPUT: stim( stim dimension[pixels or features], n images [could be ordered in time]);
%               spk( firing rate per image, n cells)
%               alpha: regularization cutoff ratio 0-1, for non-white stimuli. 
%                          0: no reg = STA, 1: fully reg = linear regression
%               delay: in frame, for stim ordered in time, STA can be computed at different time delay
% OUTPUT: sta(ndim,ncell)
%
% References:
% 1. Darragh Smyth, Ben Willmore, Gary E. Baker, Ian D. Thompson, David J. Tolhurst 
%     The Receptive-Field Organization of Simple Cells in Primary Visual Cortex of Ferrets under Natural Scene Stimulation
%     Journal of Neuroscience 1 June 2003, 23 (11) 4746-4759; 
%     https://doi.org/10.1523/JNEUROSCI.23-11-04746.2003
% 
% 2. Gidon Felsen, Jon Touryan, Feng Han, Yang Dan
%     Cortical Sensitivity to Visual Features in Natural Scenes
%     PLOS Biology 2005 
%     https://doi.org/10.1371/journal.pbio.0030342
% 
% 3. https://en.wikipedia.org/wiki/Spike-triggered_average

if nargin<3
    alpha = 0;
end

if nargin<4
    delay=0;
end

if delay>0;
    spk = [spk(delay+1:end); zeros(delay,1)];
end

dim = size(stim,1);

spk = bsxfun(@minus,  spk, mean(spk));

sta = stim*spk;

if alpha>0
    cov_stim = stim*stim'; % cov(stim'); 
    [v,d]=eig(cov_stim);
    tmp=diag(d);
    reg_rank=round(dim*alpha);
    tmp(1:end-reg_rank+1)=tmp(end-reg_rank+1);
    cov_stim_inv_partial=v*diag(1./tmp)*v';
    sta=(cov_stim_inv_partial*sta);
end

end

