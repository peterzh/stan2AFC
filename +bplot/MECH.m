function phat=MECH(w,firing_rate)
%Parameters must be column vectors, contrast must be row vector
ZL = w(:,1) + w(:,3:2:end)*firing_rate;
ZR = w(:,2) + w(:,4:2:end)*firing_rate;
phat = cat(3, exp(ZL)./(1+exp(ZL)+exp(ZR)),...
    exp(ZR)./(1+exp(ZL)+exp(ZR)),...
    1./(1+exp(ZL)+exp(ZR)) ) ;
end