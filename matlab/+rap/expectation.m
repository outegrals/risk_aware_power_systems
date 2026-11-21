function m=expectation(v)
w=[.2 .6 .2];m=reshape(sum(v.*reshape(w'*w,3,3,1),[1 2]),1,[]);
end
