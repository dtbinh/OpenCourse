function grad = computeGrad2D_neumannBCs(u)

[numRows,numCols] = size(u);

gradx = u(2:end,:) - u(1:end-1,:);
gradx = [gradx ; zeros(1,numCols)];

grady = u(:,2:end) - u(:,1:end-1);
grady = [grady,zeros(numRows,1)];

grad = cat(3,gradx,grady);

end

