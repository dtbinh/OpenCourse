function out = prox1Norm(x,t)

mx = max(abs(x) - t, 0);
out = sign(x).*mx;

end