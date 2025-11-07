function temp_Mat_sparsified = fun_far_field_sparsification(temp_Mat,threshold)

temp_Mat_sparsified = temp_Mat;

tmp = temp_Mat(:);

tmp_pos = tmp(tmp>0);
tmp_mean_pos = sum(tmp_pos)/numel(tmp_pos);


for ii=1:size(temp_Mat,1)
    temp_Mat_sparsified(ii,abs(temp_Mat(ii,:))<threshold*tmp_mean_pos) = 0;
end




% % temp_Mat_sparsified = temp_Mat;
% % 
% % for ii=1:size(temp_Mat,1)
% %     
% %     row_ii = temp_Mat(ii,:);
% %     
% %     row_ii_abs = abs(row_ii);
% %     
% %     row_ii_max = max(row_ii_abs);
% %     
% %     temp_Mat_sparsified(ii,row_ii_abs<threshold*row_ii_max) = 0;
% % end


