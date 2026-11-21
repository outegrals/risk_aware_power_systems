function report=validate_communication(out,mode)
if nargin<2,mode='full';end
root=fileparts(fileparts(mfilename('fullpath')));gold=jsondecode(fileread(fullfile(root,'output/communication/communication_validation.json')));cases=struct([]);rng(20261118,'twister');
for c=1:numel(gold.cases)
 old=gold.cases(c);n=old.n_operators;if strcmp(mode,'smoke') && n>10,continue;end
 g=rap.graph(n,old.graph);p=old.input_data.prime;b=rap.budget(g,p);k=b.iterations;assert(k==old.budget.public_iterations);
 data=old.input_data;initial=[];totals=[];
 for t=1:numel(data.runs)
  y=rap.mask(data.forecasts_encoded,data.mask_edges+1,data.runs(t).public_test_masks,p);assert(isequal(y,data.runs(t).masked_values(:)));initial(:,t)=[y;0];totals(t)=sum(y); %#ok<AGROW>
 end
 samples=unique([0 round(logspace(0,log10(k),min(k,90))) k]);trace=[];z=initial;full=k<=12000;
 if full,steps=0:k;else,steps=samples;end
 for j=steps
  if ~full
   if j==0,z=initial;elseif strcmp(g.kind,'star'),z=repmat(mean(initial),g.m,1);z(1:n,:)=z(1:n,:)+(n/(n+1))^j*(initial(1:n,:)-mean(initial(1:n,:)));elseif strcmp(g.kind,'complete'),z=repmat(mean(initial),g.m,1);else
    ev=ones(g.m,1);for s=g.shifts,ev=ev-(1-cos(2*pi*(0:g.m-1)'*s/g.m))/g.denominator;end;z=real(ifft(fft(initial).*ev.^j));
   end
  end
  if ismember(j,samples),trace(end+1,:)=[j max(abs(g.m*z-totals),[],'all')];end %#ok<AGROW>
  if full && j<k,z=rap.step(z,g);end
 end
 assert(all(mod(round(g.m*z),p)==sum(data.forecasts_encoded),'all'));
 repeats=3;times=zeros(1,repeats);executed=min(k,2000);
 for t=1:repeats,zz=initial;clock=tic;for j=1:executed,zz=rap.step(zz,g);end;times(t)=toc(clock);end
 networks=struct([]);episodes=2000;if strcmp(mode,'smoke'),episodes=100;end
 for loss=[0 .05 .2]
  for peer=[false true]
   for slot=[1 10],networks=rap.append(networks,rap.network(g,k,loss,peer,episodes,slot));end
  end
 end
 cases=rap.append(cases,struct('n',n,'graph',g.kind,'B',data.B,'delta',data.delta_MW,'budget',b,'logical_messages',2*size(g.edges,1)*k,'full_iteration_execution',full,'trace',trace,'final_error',max(abs(g.m*z-totals),[],'all'),'timed_iterations',executed,'measured_kernel_seconds',median(times),'network',networks));
 fprintf('Communication n=%d %s B=%d k=%d\n',n,g.kind,data.B,k);
end
% Event-level retry replay: receive each neighbor's frozen state before updating.
g=rap.graph(10,'circulant4');z=(1:11)';ref=z;attempts=0;
for j=1:120
 received=nan(11);pending=false(11);for e=1:size(g.edges,1),pending(g.edges(e,1),g.edges(e,2))=true;pending(g.edges(e,2),g.edges(e,1))=true;end
 while any(pending,'all'),ix=find(pending);attempts=attempts+numel(ix);success=ix(rand(numel(ix),1)>=.2);[src,~]=ind2sub([11 11],success);received(success)=z(src);pending(success)=false;end
 next=z;for s=g.shifts,for i=1:11,next(i)=next(i)+received(mod(i-1-s,11)+1,i);end;end;next=next/g.denominator;ref=rap.step(ref,g);assert(isequal(next,ref));z=next;
end
report=struct('status','passed','mode',mode,'cases',cases,'retry_replay_rounds',120,'retry_attempts',attempts,'scope','Public test masks; synthetic network times; timing uses at most 2000 actual kernel updates, never mislabeled full-run time. Network Monte Carlo uses MATLAB RNG.');rap.save_report(out,'communication',report);
f=figure('Visible','off','Color','w');hold on;for i=1:numel(cases),if cases(i).n==100 && cases(i).B==10000,semilogy(cases(i).trace(:,1),max(1e-10,cases(i).trace(:,2)),'DisplayName',cases(i).graph);end;end;set(gca,'XScale','log','YScale','log','Color','w','XColor','k','YColor','k');ylim([1e-10 1e9]);grid on;xlabel('Rounds');ylabel('Scaled error');legend('show','TextColor','k','Color','w');set(findall(f,'Type','text'),'Color','k');exportgraphics(f,fullfile(out,'communication.png'));close(f);
end
