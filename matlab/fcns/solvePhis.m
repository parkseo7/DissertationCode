function output = solvePhis(Omega, phi0, options)
% Given an initial array of phases phi0, and sync frequency, using
% step-wise iteration, locates the nearest array of phases such that the
% right-side sum of the fixed-point equation is invariant and equal to the
% frequency.

% Parameters
N = options.N;
g = options.g;
w0 = options.omega0;
gain = options.gain;
tau0 = options.tau0;

% Search options
num_iter = options.num_iter;
learning = options.learning;

% Quadratic error function (x = phi)
diff = Omega - w0;
Delta = @(phi) bsxfun(@minus, phi, phi');
tauE = @(Delta0) max(tau0+gain*sin(Delta0), 0);
E = @(phi) 0.5*sum((diff - (g/N)*sum(sin(-Omega*tauE(Delta(phi)) + Delta(phi)),2)).^2);
gradE = @(phi,i) (-g/N)*(diff - (g/N)*sum(sin(-Omega*tauE(phi - phi(i)) + phi - phi(i)))) ...
    .*sum(cos(-Omega*tauE(phi - phi(i)) + phi - phi(i)) ...
    .*(Omega*gain*cos(phi-phi(i)).*double(tauE(phi - phi(i)) > 0)-1)) ...
    + (-g/N)*sum((diff - (g/N)*sum(sin(-Omega*tauE(phi(i) - phi) + phi(i) - phi))) ...
    .*cos(-Omega*tauE(phi(i) - phi) + phi(i) - phi) ...
    .*(Omega*gain*cos(phi(i)-phi).*double(tauE(phi(i) - phi) > 0)-1));


% Iterate:
phi = phi0;
E_min = 10;
phi_min = phi0;
gradE_min = zeros(1,N);

for k=1:num_iter
    
    % Compute E:
    E_k = E(phi);
    gradE_k = zeros(1,N);
    
    for i=1:N
        gradE_k(i) = gradE(phi,i);
    end
    
    % Adjust phase (with maximum gradient)
    [MAX, ind] = max(abs(gradE_k));
    
    phi(ind) = phi(ind) - learning*E_k*phi(ind)*gradE_k(ind)/norm(gradE_k)/norm(phi);
    
    disp([num2str(k) ': ' num2str(E_k) ', ' num2str(MAX)]);
    
    % Save smallest error
    if E_k < E_min
        E_min = E_k;
        phi_min = phi;
        gradE_min = gradE_k;
    end
end

phi_c = phi_min;

% Output
output.phi = phi_c;
output.Omega = Omega;
output.err = E_min;
output.E = E;
output.gradE = gradE_min;


end

