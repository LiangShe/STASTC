function [ stc ] = STCcorr( stim, spk, alpha, delay )
%STCCORR 
% INPUT: stim(resolution, time);
%              alpha: correct ratio 0-1, for natural stimuli.
%

if nargin<4
    delay=0;
end

if delay>0;
    spk=[spk(delay+1:end); zeros(delay,1)];
end

res=size(stim,1);
DIMS=res;
ind=find(spk); 
sts=zeros(sum(spk),DIMS,'single'); % spike triggered stimuli
n_s=1;
for i=1:length(ind)
    num=spk(ind(i));
    while num>0
        sts(n_s,:)=stim(:,ind(i));
        n_s=n_s+1;
        num=num-1;
    end
end


if alpha>0 && alpha<1
    %corr
    cov_stim=cov(stim');
    [evector,evalue]=eig(cov_stim);
    reg_rank=round(res*alpha);
    tmp=diag(evalue);

    tmp(1:end-reg_rank+1)=tmp(end-reg_rank+1);
    evalue_corr=diag(tmp.^-0.5);

    sts=sts*evector*evalue_corr;

%     cov_sts=cov(sts);
    cov_sts=sts'*(sts);
    [V,D]=eig(cov_sts);
    V=evector*evalue_corr*V;
    stccorr.eigenvector=V;
    stccorr.eigenvalue=diag(D);
    stc=stccorr;
    
else
    cov_sts=cov(sts);
    [V,D]=eig(cov_sts);
    stc.eigenvector=V;
    stc.eigenvalue=diag(D);
    
end


end

