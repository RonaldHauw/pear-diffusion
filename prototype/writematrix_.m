function [out] = writematrix_(data, filename,delimiter)
%WRITEMATRIX Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(filename,'wt');
% for first row until one before last one
for row = 1:size(data,1)-1
    for col = 1:size(data, 2)-1
        fprintf(fid, '%g', data(row, col)); 
        fprintf(fid, delimiter); 
    end
    % last column
    col = size(data, 2); 
    fprintf(fid, '%g', data(row, col)); 
    fprintf(fid, '\n');
    
end
out = 1; 
% last row 
row = size(data, 1); 
for col = 1:size(data, 2)-1
    fprintf(fid, '%g', data(row, col)); 
    fprintf(fid, delimiter); 
end
col = size(data, 2); 
fprintf(fid, '%g', data(row, col)); 
fclose(fid);

end

