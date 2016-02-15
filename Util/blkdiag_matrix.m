function y = blkdiag_matrix(entries)
% Constructs a block diagonal using a 2D array of column vectors (for 1D
% block diagonal elements) or a 3d array of matrices (for 2D block diagonal
% elements).
%
% Assumes that all terms on the diagonal are of uniform size.
% Modified from the Matlab built-in function blkdiag (which takes
% non-uniform input as a sequence of arguments).
%
%                                             |A 0 .. 0|
%   Y = blkdiag_matrix([A B .. Z])  produces  |0 B .. 0|
%                                                ..
%                                             |0 0 .. Z|
%
%   Class support for inputs: 
%      float: double, single
%
%   See also DIAG, HORZCAT, VERTCAT


assert(ismember(ndims(entries),[2 3]), ...
    ['Please ensure entries is a 2D array of column vectors or a 3D array of ' ...
    'matrices. Use blkdiag for non-uniform input sizes'])

if ndims(entries) == 2
    p = 1:size(entries,1)*size(entries,2);
    m = reshape(repmat(1:size(entries,2),size(entries,1),1),1,[]);
else
%     p = repmat(reshape(1:size(entries,2)*size(entries,3),size(entries,2),size(entries,3))...
%         ,1,1,size(entries,3));
    p = repmat(reshape(1:size(entries,1)*size(entries,3),size(entries,1),size(entries,3))...
        ,1,1,size(entries,2));

    p = reshape(permute(p,[ 1 3 2 ]),size(p,1)*size(p,3),[]);
    p = p(:)';
    m = reshape(repmat(1:size(entries,2)*size(entries,3),size(entries,1),1),1,[]);
end

inds = sub2ind([p(end) m(end)],p,m);

% Get entry sizes

y = zeros(p(end),m(end));
y(inds) = entries(:);

end
