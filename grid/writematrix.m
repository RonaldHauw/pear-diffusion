function [out] = writematrix(data, filename,delimiter)
%WRITEMATRIX Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(filename,'wt');
for row = 1:size(data,1)-1
    for col = 1:size(data, 2)-1
        fprintf(fid, '%g\t', data(row, col)); 
        fprintf(fid, delimiter); 
    end
    fprintf(fid, '%g\t', data(row, col+1)); 
    fprintf(fid, '\n');
    for col = 1:size(data, 2)-1
        fprintf(fid, '%g\t', data(row+1, col)); 
        fprintf(fid, delimiter); 
    end
    fprintf(fid, '%g\t', data(row+1, col+1)); 
fclose(fid);
end
out = 1; 

end

