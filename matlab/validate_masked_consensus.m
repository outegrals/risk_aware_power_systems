function report=validate_masked_consensus(out)
root=fileparts(fileparts(mfilename('fullpath')));g=jsondecode(fileread(fullfile(root,'output/privacy/masked_consensus_validation.json')));cases=[g.correctness_cases;g.exact_rational_wraparound_case];results=struct([]);
for t=1:numel(cases)
 c=cases(t);n=c.n_operators;m=n+1;p=c.prime;y=rap.mask(c.input_forecasts_mw,c.mask_graph_edges+1,c.public_test_masks,p);assert(isequal(y,c.masked_initial_values(:)));z=[y;0];
 for k=1:c.iterations,z=[(n*z(1:n)+z(end))/m;sum(z)/m];end
 assert(all(mod(round(m*z),p)==sum(c.input_forecasts_mw)));
 if n==3
  % Exact rational recurrence with Java BigInteger (no symbolic toolbox).
  nums=arrayfun(@(v)javaObject('java.math.BigInteger',sprintf('%.0f',v)),[y;0],'UniformOutput',false);den=javaObject('java.math.BigInteger','1');N=javaObject('java.math.BigInteger',num2str(n));M=javaObject('java.math.BigInteger',num2str(m));four=javaObject('java.math.BigInteger','4');
  for k=1:c.iterations
   next=nums;total=javaObject('java.math.BigInteger','0');for i=1:m,total=total.add(nums{i});end
   for i=1:n,next{i}=nums{i}.multiply(N).add(nums{m});end
   next{m}=total;nums=next;den=den.multiply(M);
  end
  sumY=javaObject('java.math.BigInteger',sprintf('%.0f',sum(y)));
  for i=1:m,diff=nums{i}.multiply(M).subtract(den.multiply(sumY)).abs();assert(diff.multiply(four).compareTo(den)<0);end
 end
 results=rap.append(results,struct('n',n,'iterations',c.iterations,'max_error',max(abs(m*z-sum(y))),'exact_integer_check',n==3));
end
tri=[1 2;1 3;2 3];a=[1 1 2];b=[1 2 1];
assert(isequal(views(a,tri,[]),views(b,tri,[])));assert(isequal(views(a,tri,1),views(b,tri,1)));
assert(~isequal(views(a,[1 2],[]),views(b,[1 2],[])));assert(~isequal(views(a,[1 2;1 3],1),views(b,[1 2;1 3],1)));
assert(isequal(([120 100 80]*.75)/.75,[120 100 80]));
report=struct('status','passed','cases',results,'exhaustive_triangle_masks_per_input',343,'unmasked_attack','initial and first star-leaf states reveal forecasts','scope','Finite passive-adversary checks, public test masks; not deployed cryptography');rap.save_report(out,'masked_consensus',report);
end
function V=views(x,e,colluders)
r=rap.product_states(7,size(e,1));known=any(ismember(e,colluders),2);V=zeros(size(r,1),numel(colluders)+sum(known)+3);
for j=1:size(r,1),V(j,:)=[x(colluders),r(j,known),rap.mask(x,e,r(j,:),7)'];end
V=sortrows(V);
end
