function out = proxHuber(x,mu,t)

% We compute prox_{tf}(x) where f is the Huber penalty with parameter mu.

% out = (t*prox1Norm(x,t) + mu*x)/(t + mu); % Shouldn't the step size in prox1Norm be t + mu ???

% out = (t*prox1Norm(x,t + mu) + mu*x)/(t + mu); % Shouldn't the step size in prox1Norm be t + mu ???

out = (t.*prox1Norm(x,t + mu) + mu*x)./(t + mu); % Shouldn't the step size in prox1Norm be t + mu ???

end