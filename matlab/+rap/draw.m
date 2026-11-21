function v=draw(p,key,indices)
n=numel(indices)-2;ix=sub2ind([3 3 n],repmat(indices(1),1,n),indices(2:n+1),1:n);a=p.(key);v=reshape(a(ix),1,[]);
end
