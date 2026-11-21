function value=upward_power(x,k,mc)
value=javaObject('java.math.BigDecimal','1');
while k>0
 if mod(k,2),value=value.multiply(x,mc);end
 k=floor(k/2);if k>0,x=x.multiply(x,mc);end
end
end
