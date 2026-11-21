function a=assets(f)
a.size=f.size(:)';a.phase=f.phase(:)';a.capacity=4*a.size;a.power=.65*a.size;a.export_limit=1.8*a.size;
end
