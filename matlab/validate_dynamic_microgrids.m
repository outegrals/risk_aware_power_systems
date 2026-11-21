function report=validate_dynamic_microgrids(out,mode)
if nargin<2,mode='full';end
root=fileparts(fileparts(mfilename('fullpath')));fixtures=jsondecode(fileread(fullfile(root,'matlab/fixtures/microgrid_inputs.json')));
runs=struct([]);trace=[];r=0;start=tic;
for n=[10 50 100]
 for scale=[0 .5 1 2]
  reps=0:29;if strcmp(mode,'smoke'),reps=0;end
  for rep=reps
   f=fixtures([fixtures.n]==n & [fixtures.rep]==rep);a=rap.assets(f);energy=.5*a.capacity;expected=zeros(3,1);cash=expected;shortfall=expected;esf=expected;schedule=expected;hits=expected;
   margins=inf(1,4);maxres=0;balanceerr=0;socerr=0;imports=0;spill=0;throughput=0;soclo=1;sochi=0;timings=zeros(192,1);
   for step=0:191
    p=rap.physics(a,energy,step,scale);[mu,u,C,load]=rap.inputs(p,5/n);t=tic;[q,k,res]=rap.strategies(mu,u);timings(step+1)=toc(t);maxres=max(maxres,res);Q=sum(q,2);
    expected=expected+.25*(sum(C)+2*(5/n)*Q.*(sum(mu)-Q));schedule=schedule+.25*Q;hits=hits+sum(abs(q-u)<1e-10,2);
    for j=1:3,esf(j)=esf(j)+.25*sum(rap.expectation(max(reshape(q(j,:),1,1,n)-p.export,0)));end
    bal=p.generation+p.discharge+p.import-p.load-p.charge-p.export-p.spill;balanceerr=max(balanceerr,max(abs(bal),[],'all'));
    ev=max([0;reshape(.1*reshape(a.capacity,1,1,n)-p.next_energy,[],1);reshape(p.next_energy-.9*reshape(a.capacity,1,1,n),[],1)]);socerr=max(socerr,ev);
    assert(max(p.charge.*p.discharge,[],'all')<1e-10);assert(max(p.charge+p.discharge-reshape(a.power,1,1,n),[],'all')<1e-10);
    signed=p.export-p.import;lo=min(p.external_cond_mean-.02*scale*p.external_base-reshape(sum(max(signed,[],2),3),3,1));hi=max(p.external_cond_mean+.02*scale*p.external_base-reshape(sum(min(signed,[],2),3),3,1));cap=5*sum(a.size);
    margins=min(margins,[load-sum(u),cap-load,lo,cap-hi]);assert(min(margins)>0 && ev<1e-10);
    ix=double(f.draws_zero_based(step+1,:))+1;x=rap.draw(p,'export',ix);im=rap.draw(p,'import',ix);ch=rap.draw(p,'charge',ix);dis=rap.draw(p,'discharge',ix);sp=rap.draw(p,'spill',ix);
    z=[-1 0 1];ext=p.external_cond_mean(ix(1))+.02*scale*p.external_base*z(ix(end));rt=20+2*(5/n)*(ext+sum(im)-sum(x));da=20+2*(5/n)*(load-Q);
    cash=cash+.25*(da.*Q+rt*(sum(x-im)-Q)-5*sum(ch+dis));shortfall=shortfall+.25*sum(max(q-x,0),2);
    imports=imports+.25*sum(im);spill=spill+.25*sum(sp);throughput=throughput+.25*sum(ch+dis);energy=rap.draw(p,'next_energy',ix);soclo=min(soclo,min(energy./a.capacity));sochi=max(sochi,max(energy./a.capacity));
    if n==50 && scale==1 && rep==0,trace(end+1,:)=[(step+1)*.25 mean(energy./a.capacity) min(energy./a.capacity) max(energy./a.capacity) mean(x) mean(im) mean(mu) mean(u) mean(q,2)' mean(energy)];end %#ok<AGROW>
   end
   r=r+1;runs=rap.append(runs,struct('n',n,'uncertainty_scale',scale,'replicate',rep,'expected_cash',expected','realized_cash',cash','scheduled_MWh',schedule','shortfall_MWh',shortfall','expected_shortfall_MWh',esf','cap_hit_fraction',(hits/(n*192))','imports_MWh',imports,'curtailment_MWh',spill,'throughput_MWh',throughput,'soc_range',[soclo sochi],'initial_energy_MWh',sum(.5*a.capacity),'terminal_energy_MWh',sum(energy),'dispatch_margins_MW',margins,'max_br_residual',maxres,'max_balance_residual',balanceerr,'max_soc_violation',socerr,'solver_median_seconds',median(timings)));
  end
  fprintf('Microgrids n=%d scale=%g complete\n',n,scale);
 end
end
% Compare all published trajectory quantities on identical random inputs.
golden=jsondecode(fileread(fullfile(root,'output/microgrid/results.json')));maxerr=0;
fields={'expected_cash','realized_cash','scheduled_MWh','shortfall_MWh','expected_shortfall_MWh','imports_MWh','curtailment_MWh','throughput_MWh','terminal_energy_MWh'};
for i=1:numel(runs)
 g=golden.runs([golden.runs.n]==runs(i).n & [golden.runs.uncertainty_scale]==runs(i).uncertainty_scale & [golden.runs.replicate]==runs(i).replicate);
 for j=1:numel(fields),v=runs(i).(fields{j});w=g.(fields{j});err=max(abs(v(:)-w(:)));assert(err<1e-7*max(1,max(abs(w(:)))));maxerr=max(maxerr,err);end
end
summary=struct([]);j=0;
for n=[10 50 100]
 for s=[0 .5 1 2]
  a=runs([runs.n]==n & [runs.uncertainty_scale]==s);cash=vertcat(a.realized_cash)/n;gain=cash(:,3)-cash(:,1);ci=[NaN NaN];if numel(a)==30,margin=2.04523*std(gain)/sqrt(30);ci=mean(gain)+[-margin margin];end
  j=j+1;summary=rap.append(summary,struct('n',n,'scale',s,'cash_per_operator',mean(cash,1),'paired_95pct_t_interval',ci,'shortfall_kWh_per_operator',1000*mean(vertcat(a.shortfall_MWh),1)/n));
 end
end
report=struct('status','passed','mode',mode,'rng','Recorded Python random inputs; all dynamics and calculations run natively in MATLAB','max_absolute_python_parity_error',maxerr,'intervals',192*numel(runs),'elapsed_seconds',toc(start),'summary',summary,'runs',runs);
rap.save_report(out,'dynamic_microgrids',report);
T=array2table(trace,'VariableNames',{'hour_end','mean_soc','min_soc','max_soc','mean_export_MW','mean_import_MW','mean_forecast_MW','mean_cap_MW','RN_mean_bid_MW','ERS_mean_bid_MW','retuned_mean_bid_MW','mean_next_energy_MWh'});writetable(T,fullfile(out,'representative_trajectory.csv'));
end
