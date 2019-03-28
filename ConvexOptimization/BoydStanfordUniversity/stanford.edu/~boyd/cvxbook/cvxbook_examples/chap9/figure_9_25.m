% Generates figure 9.25, Boyd & Vandenberghe, Convex Optimization
%
% Number of Newton iterations required to minimize self-concordant
% functions.

warning off

% analytic center of  Ax <= b, where b>0 
randn('seed',0);
rand('seed',0);
mvalues = [100, 1000, 1000];
nvalues = [ 50,  500,   50];
noproblems = 5;       % number of problems for each dimension 
nostartingpts = 10;   % number of starting points for each problem

niters_1 = []; initvals_1 = [];
niters_2 = []; initvals_2 = [];
niters_3 = []; initvals_3 = [];

alpha = 0.1;
beta = 0.8;
for probsize = 1: length(mvalues)
    m = mvalues(probsize);  
    n = nvalues(probsize);
    prob = 0; 
    while prob < noproblems

        A = randn(m,n);   
        b = rand(m,1);

        % Find exact minimum.

        x = zeros(n,1);
        maxiters = 20;
        unbounded = 1;  
        for iters=1:maxiters
            y = b-A*x;
            g = A'*(1./y);
            H = A'*diag(1./y.^2)*A;
            v = -H\g;
            lambda = sqrt(-g'*v);
            if (lambda < 1e-8), unbounded = 0;  break; end;
            dy = -A*v;
            s = 1; 
            ysc = dy./y;
            while (min(1+s*ysc) < 0), s = beta*s; end;
            f = -sum(log(y));
            fnew  = f - sum(log(1+s*ysc));
            while (fnew > f - alpha*s*lambda^2),  
                s = beta*s; 
                fnew  = f - sum(log(1+s*ysc));
            end;
            x = x+s*v;
        end;
        if (unbounded == 1) 
            disp('Unbounded.'), continue;   % generate another problem
        else 
            prob = prob+1; 
        end;
        xmin = x;
        fmin = -sum(log(b-A*xmin));
        y = b-A*xmin;  

        for k = 1:nostartingpts

            % Generate a random direction v and compute
            %
            %     smax = sup{s | xmin + sv > 0} 
            %     fmax = f(xmin + (1-1e-6)*smax*v).
            %
            % Take a starting point x = xmin + s*v such that 
            %
            %     f(x) = fmin + (k/nostartingpoints)*(fmax - fmin).

            v = randn(n,1);  
            dy = -A*v;
            inds = find(dy<0); 
            smax = min( -y(inds)./dy(inds) );
            maxinitdiff = -sum( log(1+(1-1e-6)*smax*dy./y) );
            initdiff = (k/nostartingpts) * maxinitdiff;
            l = 0;  u = smax;
            s = exp((log(l+1)+log(u+1))/2)-1;
            f0 = -sum(log(1+s*dy./y));
            while ((abs(f0 - initdiff) > 1e-3) & (u-l > 1e-10))
                if f0 > initdiff, 
                    u=s; 
                else 
                    l=s;  
                end;
                s = exp((log(l+1)+log(u+1))/2)-1;
                f0 = -sum(log(1+s*dy./y));
            end;
            x = xmin + s*v;
            if (probsize==1) 
                initvals_1 = [initvals_1, -sum(log(b-A*x)) - fmin];
            elseif (probsize==2) 
                initvals_2 = [initvals_2, -sum(log(b-A*x)) - fmin];
            else  
                initvals_3 = [initvals_3, -sum(log(b-A*x)) - fmin];
            end; 

            % Run Newton's method.
            maxiters = 100;
            for iters=1:maxiters
                y = b-A*x;
                g = A'*(1./y);
                H = A'*diag(1./y.^2)*A;
                v = -H\g;
                if (g'*v > 0), iters = Inf; break; end; 
                lambda = sqrt(-g'*v);
                if (lambda < 1e-5), break; end;
                dy = -A*v;
                s = 1; 
                ysc = dy./y;
                while (min(1+s*ysc) < 0), s = beta*s; end;
                f = -sum(log(y));
                fnew  = f - sum(log(1+s*ysc));
                while (fnew > f - alpha*s*lambda^2),  
                    s = beta*s; 
                    fnew  = f - sum(log(1+s*ysc));
                end;
                x = x+s*v;
            end;
            if (iters == maxiters), iters = Inf; end;
            disp(['Problem size ', int2str(probsize), ', problem ',...
                int2str(prob), ', starting point ', int2str(k), ...
                '.  Number of Iterations: ', int2str(iters), '.']); 
            if (probsize==1) 
                niters_1 = [niters_1, iters];
            elseif (probsize==2) 
                niters_2 = [niters_2, iters];
            else  
                niters_3 = [niters_3, iters]; 
            end;
        end;
    end;
end;

d = plot(initvals_1, niters_1, 'o', initvals_2, niters_2, 's', ...
     initvals_3, niters_3, 'd');
xlabel('f(x0)-pstar');
ylabel('iterations');
axis([0 35 0 26]);  
