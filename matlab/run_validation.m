function summary=run_validation(mode)
% RUN_VALIDATION('full') executes all studies. 'smoke' runs one dynamic replicate
% per setting and small communication graphs. No manuscript is modified.
if nargin<1,mode='full';end
assert(ismember(mode,{'full','smoke'}));root=fileparts(fileparts(mfilename('fullpath')));out=fullfile(root,'output/matlab',mode);if ~isfolder(out),mkdir(out);end
start=tic;summary=struct('mode',mode,'matlab',version,'computer',computer,'status','running');
validate_risk_and_caps(out);validate_masked_consensus(out);validate_incentives(out);
validate_dynamic_microgrids(out,mode);audit_dynamic_microgrids(out);validate_dispatch_feasibility(out);render_microgrid_evidence(out);validate_communication(out,mode);
summary.status='passed';summary.elapsed_seconds=toc(start);rap.save_report(out,'validation_summary',summary);fprintf('All MATLAB validations passed: %s\n',out);
end
