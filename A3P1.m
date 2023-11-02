clear 
close all

cd HW3_DataF23
W = importdata('WordDocF23.dat'); % 12x10 (12 words and 10 docs)
cd ..

% SVD (singular value decomposition)
[U,S,V] = svd(W,"econ");


% % Perform PCA to W:
% Covariance_matrix = cov(W); % Get the covariance matrix
% [eigenvector_matrix,eigenvalue_matrix] = eig(Covariance_matrix);
% eigenvalues = diag(eigenvalue_matrix);
% max_eigenvalue = max(eigenvalues);

% k=1 approximation:
k = 1;
W1 = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
error1 = (norm(W1 - W, 'fro'))^2;

% k=3 approximation:
k = 3;
W3 = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
error3 = (norm(W3 - W, 'fro'))^2;

% k=5 approximation:
k = 5;
W5 = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
error5 = (norm(W5 - W, 'fro'))^2;

% The most similar document pair: (rank-3 approximation)
Doc_Inner_similarity_matrix = W3'* W3;
% Extract upper triangle without diagonal
upper_triangle_doc = triu(Doc_Inner_similarity_matrix, 1); % Extracts upper triangle excluding the diagonal
% Find the maximum value and its indices
[max_value_doc, ind] = max(upper_triangle_doc, [], 'all', 'linear');
% Convert linear index to row and column indices
[row_doc, col_doc] = ind2sub(size(Doc_Inner_similarity_matrix), ind);
% % Display the max value and its indices
% disp(max_value);
% disp([row, col]);


% The most similar words pair: (rank-3 approximation)
Word_Inner_similarity_matrix = W3* W3';
% Extract upper triangle without diagonal
upper_triangle_word = triu(Word_Inner_similarity_matrix, 1); % Extracts upper triangle excluding the diagonal
% Find the maximum value and its indices
[max_value_word, ind] = max(upper_triangle_word, [], 'all', 'linear');
% Convert linear index to row and column indices
[row_word, col_word] = ind2sub(size(Word_Inner_similarity_matrix), ind);
% % Display the max value and its indices
% disp(max_value);
% disp([row, col]);



% Finding cos similarity:
% Doc:
V3_bar = S(1:3,1:3)*(V(:,1:3))';
% V3_bar'*V3_bar will be the same as doc inner product similarity
Doc_cos_similarity_matrix = zeros(10,10);
for i = 1:10
    for j = 1:10
        Vi = V3_bar(:,i);
        Vj = V3_bar(:,j);
        Doc_cos_similarity_matrix(i,j) = (Vi'*Vj)/(norm(Vi)*norm(Vj));
    end
end
upper_triangle_cos_doc = triu(Doc_cos_similarity_matrix, 1); % Extracts upper triangle excluding the diagonal
[max_cos_doc, ind] = max(upper_triangle_cos_doc, [], 'all', 'linear');
[row_cos_doc, col_cos_doc] = ind2sub(size(Doc_cos_similarity_matrix), ind);

% Word:
U3_bar = U(:,1:3)*S(1:3,:);
U3_bar*U3_bar'
% U3_bar*U3_bar' will be the same as word inner product similarity
Word_cos_similarity_matrix = zeros(12,12);
for i = 1:12
    for j = 1:12
        Ui = U3_bar(i,:);
        Uj = U3_bar(j,:);
        Word_cos_similarity_matrix(i,j) = (Ui*Uj')/(norm(Ui)*norm(Uj));
    end
end
upper_triangle_cos_word = triu(Word_cos_similarity_matrix, 1); % Extracts upper triangle excluding the diagonal
[max_cos_word, ind] = max(upper_triangle_cos_word, [], 'all', 'linear');
[row_cos_word, col_cos_word] = ind2sub(size(Doc_cos_similarity_matrix), ind);
