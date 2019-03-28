% Matlab script for invoking sdpsol from within Matlab
% Require Matlab string variable:
%   SDPSOL_FILENAME -- contains sdpsol source filename
%
% last modified: 05/13/96

save .sdpsol_tmp.mat
eval(['!sdpsol -I .sdpsol_tmp.mat -m .sdpsol_out.mat ' SDPSOL_FILENAME])
load .sdpsol_out.mat
!'rm' .sdpsol_tmp.mat
!'rm' .sdpsol_out.mat
