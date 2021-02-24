function noiseFun = noiseMollifier(t1, t2)
% Returns the noise function for finite-time stimulation on [t1, t2]
% Ensures the noise is smooth.

% Smooth mollifier
tp = (t1 + t2) / 2;
COEFF = 1 / exp(-(tp-t1).^-2).*exp(-(tp-t2).^-2);
MOL = @(t) COEFF*exp(-(t-t1).^-2).*exp(-(t-t2).^-2);

OnInterval = @(t) double(t > t1).*double(t < t2);

noiseFun = @(t) OnInterval(t).*MOL(t);

end


