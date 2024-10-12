dimension = 2;
function [val, df] = fdf(x)
  dim = size(x)(1);
  vec = zeros(1,dim);
  cost = 0;
  for i = 1:dim-1,
      cost = cost + (1 - x(i))^2 + 100*(x(i+1) - x(i)^2)^2;
      vec(i) = -2*(1 - x(i)) - 400*x(i)*(x(i+1) - x(i)^2);
  end
  vec(dim) = 200*(x(dim) - x(dim-1)^2);

  val = cost;
  df = vec;
endfunction

// fonction h retournant la hessienne  de f
function [dff] = Hv(x,v)
  dim = size(v)(1);
  G = zeros(dim,dim);
  for i = 1:dim-1,
      G(i,i) = 2 - 400*(x(i+1) - 3*x(i)^2);
      G(i,i+1) = -400*x(i);
  end
  G(dim,dim) = 200;
  G(dim,dim-1) = -400*x(dim-1);
  
  dff = v'*G;
endfunction

verbose = 1;
MaxIter = 20;
exec ("GC_TR.sci",-1);
c = zeros(1,dimension);
Delta = 0.9;
x0 = rand(1,dimension);
[dNTR,dNQ,ngc] = GC_TR(c, Hv,x0, 1e-8, MaxIter, Delta,verbose);

