function fxn_plot_cal(input1, input2, c_range)
%%

for i = 1:cell_file1_num
    for ii = 1:cell_file2_num
        jAB_inter{i,ii} = intersect(cell_file1_input{i,1},cell_file2_input{ii,1});
        jAB_unique{i,ii} = unique([cell_file1_input{i,1};cell_file2_input{ii,1}]);
%         jAB_jaccard{i,ii} = numel(jAB_inter{i,1})/numel(jAB_unique{ii,1});
    end
end

for i = 1:cell_file1_num
    for ii = 1:cell_file2_num
          jAB_jaccard{i,ii} = numel(jAB_inter{i,ii})/numel(jAB_unique{i,ii});
    end
end

jAB_jaccard_mat = cell2mat(jAB_jaccard);
zero_fill = find(jAB_jaccard_mat ==1);
jAB_jaccard_mat_zero_fill = jAB_jaccard_mat;
jAB_jaccard_mat_zero_fill(zero_fill) = 0;

%%
end