function b=budget(g,p)
% Directed-upward decimal arithmetic using Java BigDecimal; no Python required.
mc=javaObject('java.math.MathContext',int32(70),java_get('java.math.RoundingMode','CEILING'));
D=@(x)javaObject('java.math.BigDecimal',sprintf('%.0f',x));one=D(1);m=g.m;C=m*ceil(sqrt(m))*(p-1);rho=D(g.rho_num).divide(D(g.rho_den),mc);
if g.rho_num==0,k=1;else,k=max(1,ceil(log(.20/C)/log(g.rho_upper)));end
ideal=D(C).multiply(rap.upward_power(rho,k,mc),mc);quarter=javaObject('java.math.BigDecimal','0.25');
while ideal.compareTo(quarter)>=0,k=k+1;ideal=D(C).multiply(rap.upward_power(rho,k,mc),mc);end
u=one.divide(D(2^53),mc);tu=D(m+5).multiply(u,mc);gamma=tu.divide(one.subtract(tu),mc);
roundoff=D(m*(p-1)).multiply(rap.upward_power(one.add(gamma,mc),k,mc).multiply(one.add(u,mc),mc).subtract(one),mc);
total=ideal.add(roundoff,mc);assert(total.compareTo(javaObject('java.math.BigDecimal','0.5'))<0);
b=struct('iterations',k,'ideal_bound_upward',char(ideal.toString()),'roundoff_bound_upward',char(roundoff.toString()),'total_bound_upward',char(total.toString()));
end
function v=java_get(cls,field)
v=javaMethod('valueOf',cls,field);
end
