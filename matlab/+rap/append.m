function list=append(list,item)
if isempty(list),list=item;else,list(end+1)=item;end
end
