clear;clc;

x1 = -3:0.2:3;
x2 = -3:0.2:3;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];


y = zeros(size(X1));

for mu_val=[0,-1.2,2]
    mu = [mu_val mu_val];
    Sigma = [0.25 0.3; 0.3 1];

    y1 = mvnpdf(X,mu,Sigma);
    y1 = reshape(y1,length(x2),length(x1));

    y = y + y1;
end
% 
% mu = [2 2];
% 
% y2 = mvnpdf(X,mu,Sigma);
% y2 = reshape(y2,length(x2),length(x1));
% 
% y = y1 + y2;

surf(x1,x2,y)
caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
axis([-3 3 -3 3 0 1])
xlabel('x1')
ylabel('x2')
zlabel('Probability Density')
