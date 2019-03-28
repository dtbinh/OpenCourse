function [z_loc L_loc] = sens_sel_loc(A, zhat)
% Local optimization method for sensor selection as described in the paper
% Sensor Selection via Convex Optimization
% www.stanford.edu/~boyd/papers/sensor_selection.html
%
% May 2008 Siddharth Joshi & Stephen Boyd

[m n] = size(A);
%k = sum(zhat);
zn = zhat;
iteration = 0; 
N_loc = 10*m^3/n^2; 
swapstaken = 0;

hatSig = inv(A'*diag(zn)*A);

while(iteration < N_loc)
    flag = 0;
    S = find(zn == 1); nS = find(zn == 0);   
    
    for out = S' %matlab does not allow to say for out = S
        for in = nS' 
            iteration = iteration + 1;
            v1 = [A(out,:); A(in, :)];
            v2 = [-A(out,:); A(in, :)]';
            S_2by2 = eye(2) + v1*hatSig*v2;
            volchangefactor = det(S_2by2);
            if(volchangefactor > 1)
                %fprintf('Value= %f\t Swap: OUT %d\tIN %d\n', val, out, in);
                swapstaken = swapstaken + 1;
                hatSig = hatSig - hatSig*v2*inv(S_2by2)*v1*hatSig;
                zn(out) = 0; zn(in) = 1;
                flag = 1;
                break;
            end
        end
        if( flag == 1), break; end
    end
    if(flag==0), break; end
end

fprintf('\nLocal optimization:\n');
if(flag==1)
   fprintf('Maximum Iteration = %d reached. Terminated.\n', N_loc);
end
fprintf('Swaps checked = %d, swaps taken = %d, [k(m-k) = %d]\n', iteration, swapstaken, length(S)*length(nS));
z_loc = zn;
L_loc = log(det(A'*diag(z_loc)*A));

end
