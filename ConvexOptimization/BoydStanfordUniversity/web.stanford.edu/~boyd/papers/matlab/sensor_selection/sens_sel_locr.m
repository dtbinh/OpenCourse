function [z_loc L_loc] = sens_sel_locr(A, k, zast, threshold)
% Local optimization method for sensor selection as described in the paper
% Sensor Selection via Convex Optimization
% www.stanford.edu/~boyd/papers/sensor_selection.html
%
% Nov 2007 Siddharth Joshi & Stephen Boyd

if(nargin < 4); threshold = 0.5; end

[m n] = size(A);
N_loc = 10*m^3/n^2;

[z_s idx_s] = sort(zast); 
thres = z_s(m-k); zhat=(zast>thres); zhat_s = (z_s > thres);
toswap01 = (abs(z_s - 0.5) <= threshold);
S = idx_s(toswap01 & zhat_s); 
nS = flipud(idx_s(toswap01 & ~zhat_s)); 

fprintf('\nLocal optimization (ordering and thresholding):\n');
fprintf('Number of sensors to chosen to swap: %d ', sum(toswap01));
fprintf('selected: %d, not selected %d\n', length(S), length(nS));

iteration = 0; 
swapstaken = 0;

hatSig = inv(A'*diag(zhat)*A);

while(iteration < N_loc)
    flag = 0;
    
    for out = S' %matlab does not allow to say for out = S
        for in = nS' 
            iteration = iteration + 1;
            v1 = [A(out,:); A(in, :)];
            v2 = [-A(out,:); A(in, :)]';
            S_2by2 = eye(2) + v1*hatSig*v2;
            volchangefactor = det(S_2by2);
            if(volchangefactor > 1)
                %fprintf('Value= %f\t Swap: OUT %d\tIN %d\n', volchangefactor, out, in);
                swapstaken = swapstaken + 1;
                hatSig = hatSig - hatSig*v2*inv(S_2by2)*v1*hatSig;
                zhat_s(idx_s == out) = 0; zhat_s(idx_s == in) = 1;
                S = idx_s(toswap01 & zhat_s); 
                nS = flipud(idx_s(toswap01 & ~zhat_s)); 
                flag = 1;
                break;
            end
        end
        if( flag == 1), break; end
    end
    if(flag==0), break; end
end


if(flag==1)
   fprintf('Maximum Iteration = %d reached. Terminated.\n', N_loc);
end
fprintf('Swaps checked = %d, swaps taken = %d, len(S)*len(ns) = %d]\n', iteration, swapstaken, length(S)*length(nS));
z_loc = zeros(m,1); z_loc(idx_s(zhat_s)) = 1;
L_loc = log(det(A'*diag(z_loc)*A));

end
