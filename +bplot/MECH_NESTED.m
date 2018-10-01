function phat=MECH_NESTED(b_GO,b_LR,w_GO,w_LR,firing_rate)
%Parameters must be column vectors, contrast must be row vector


ZGO = b_GO + w_GO(:,1).*firing_rate(1,:) + w_GO(:,2).*firing_rate(2,:) + w_GO(:,3).*firing_rate(3,:) + w_GO(:,4).*firing_rate(4,:);
ZLR = b_LR + w_LR(:,1).*firing_rate(1,:) + w_LR(:,2).*firing_rate(2,:) + w_LR(:,3).*firing_rate(3,:) + w_LR(:,4).*firing_rate(4,:);

pGO = exp(ZGO)./(1+exp(ZGO));
pLGO = exp(ZLR)./(1+exp(ZLR));
pRGO = 1-pLGO;

phat = cat(3, pLGO.*pGO,...
    pRGO.*pGO,...
    1-pGO ) ;
end