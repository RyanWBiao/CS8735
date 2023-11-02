close all
clear

cd HW3_DataF23
Distances = importdata('MOCityDistF23.dat'); % 12x10 (12 words and 10 docs)
cd ..

% step1 Compute the distance matrix:
square_distance_matrix = zeros(9,9);
for i = 1:9
    for j=1:9
        square_distance_matrix(i,j) = Distances(i,j)*Distances(i,j);
    end
end

% step2: double centering on the distance matrix
J = eye(9) - (1/9)*ones(9,9);
B = -0.5*J*square_distance_matrix*J;

% step3: eigendecomposition on the double_centered distance matrix:
[V, D] = eig(B);
eigenvalues = diag(D);
[sorted_eigenvalues, idx] = sort(eigenvalues, 'descend');
sorted_eigenvectors = V(:, idx);

% step4: Compute the eigen components:
for k=1:8
    if sorted_eigenvalues(k)>0 && sorted_eigenvalues(k+1)<0
        break
    end
end

positive_eigenvalue_matrix = zeros(k,k);
for i = 1:k
    positive_eigenvalue_matrix(i,i) = sorted_eigenvalues(i);
end

positive_eigenvector_matrix = sorted_eigenvectors(:,1:k);

% step5: Form the coordinate matrix for the data samples:
% perserve only 2 eignevalues and eigenvetors:
eigenvalues_2 = positive_eigenvalue_matrix(1:2,1:2);
eigenvectors_2 = positive_eigenvector_matrix(:,1:2);

eigenvectors_2(:,1) = -1*eigenvectors_2(:,1); %swap direction

coordinates = zeros(9,2);
for i = 1:2
    coordinates(:,i) = sqrt(eigenvalues_2(i,i))*eigenvectors_2(:,i);
end

% scatter plot:
scatter(coordinates(:,1), coordinates(:,2));

labels = {'Branson', 'Cape Girardeau', 'Columbia', 'Jefferson City', 'Kansas City','Rolla','St. Louis','Springfield','St. Jeseph'};
for i = 1:length(coordinates)
    text(coordinates(i,1), coordinates(i,2), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Add title and labels
title('MO cities');
xlabel('X-axis');
ylabel('Y-axis');

% Save the plot as a PDF file
saveas(gcf, 'MO cities.pdf');




