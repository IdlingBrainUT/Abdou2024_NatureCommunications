function [dot_matrix] = fxn_dot_product(basis_Y, basis_X, color_range)
%% comment
% matrix shouled consist (cells x pattern) matrix.
% the number of cells should take same number.

% color_range = [0 1];
%% Calculatin dot product
% C = dot(A,B,dim)
% inner_product a*b = |a|*|b|*cos
% cos = inner_product  / |a|*|b|

% basis_Y = data_table{3,10};  % {x,10} = thresholding;
% basis_X = data_table{10,10}; % {x,10} = thresholding;

dot_matrix = zeros(size(basis_Y,2),size(basis_X,2)); % initialize map
for i_X = 1:size(basis_X,2)
for i_Y = 1:size(basis_Y,2) 
inner_p = dot(basis_Y(:,i_Y),basis_X(:,i_X));
normA = norm(basis_Y(:,i_Y));
normB = norm(basis_X(:,i_X));
cos = inner_p/(normA*normB);
    dot_matrix(i_Y,i_X) = cos;
end
end
%% Figure
% figure('Position', [100 200 200 150]);
imagesc(dot_matrix);
% title('Jaccard index')
% xlabel(['Session-', file2_name]); ylabel(['Session-', file1_name]); % colorbar;
clim([color_range]); colormap('parula');
set(gca, 'FontSize',7, 'FontName','Arial');
%%
end
%%