function p=physics(a,energy,step,scale)
% Common x local x operator, same support and bid-independent policy as reference.
hour=mod(step*.25,24);sz=reshape(a.size,1,1,[]);phase=reshape(a.phase,1,1,[]);energy=reshape(energy,1,1,[]);
wind=1.8*sz.*(.62+.10*sin(2*pi*hour/24+phase));solar=2.4*sz*max(0,sin(pi*(hour-6)/12));
base=.65*sz.*(.85+.08*cos(2*pi*hour/24+phase)+.18*exp(-((hour-19)/3)^2));
c=reshape([-1 0 1],3,1,1);l=reshape([-1 0 1],1,3,1);
p.generation=(wind+solar).*(1+scale*(.18*c+.12*l));p.load=base.*(1+scale*(.04*c+.04*l));
surplus=p.generation-p.load;ref=sz*(.60+.10*sin(2*pi*hour/24));
mc=min(.65*sz,(.9*4*sz-energy)/(.95*.25));md=min(.65*sz,.95*(energy-.1*4*sz)/.25);
p.charge=min(max(surplus-ref,0),max(0,mc));p.discharge=min(max(ref-surplus,0),max(0,md));net=surplus+p.discharge-p.charge;
p.export=min(max(net,0),1.8*sz);p.import=max(-net,0);p.spill=max(net-1.8*sz,0);p.next_energy=energy+.25*(.95*p.charge-p.discharge/.95);
p.external_base=2.7*sum(a.size)*(1+.10*cos(2*pi*(hour-19)/24));p.external_cond_mean=p.external_base*(1+.06*scale*[-1;0;1]);
end
