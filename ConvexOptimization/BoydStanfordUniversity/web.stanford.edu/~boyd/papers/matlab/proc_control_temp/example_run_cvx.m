% EXAMPLE_RUN_CVX  solves the processor speed problem using CVX
%                  "Processor speed control with thermal constraints"
%                  by Mutapcic, Boyd, Murali, Atienza, De Micheli, Gupta
%
% CVX is freely available at www.stanford.edu/~boyd/cvx

% create problem data for the example in the paper
create_example;

%
% CVX solution
%
cvx_begin
  variable s_opt(n)    % speed variable for each processor
  maximize sum(s_opt)  % maximize total processing power (throughput)
  subject to
    % maximum temperature constraint
    G*pow_pos(s_opt,3) + Tother + Tamb <= Tmax;
    % processor speed bounds
    smin <= s_opt; s_opt <= smax;
cvx_end
optval = cvx_optval;
fprintf(1,'Optimal throughput is %3.4f\n',optval);
