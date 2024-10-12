function [dN, dNQ, niter, L_f] = GC_Hv(c, Hv, x, eps, MaxIter, verbose, delt0)
// Entrée:
//    c(1,n), Q(n,n) paramètres de la fonction quadratique à minimiser
//        q(delt) = 0.5*delt'*Q*delt + c*delt
//        nabla q(delt) = delt'*Q + c
//    eps tolérance d'arrêt
//    MaxIter  nombre maximum d'itérations
//    verbose  imprime les itérations (>0) ou non (<=0)
//    delt0(n,1)  point de départ des itérations
//
// Sortie:
//    dN: la solution  :dN*Q+b=0
//    niter: nombre total d'itérations
//    L_f: liste des valeurs de fonctions au fil des itérations
  
  // init
  niter = 0; 
  L_f = [];
  [bidon,dim] = size(c)
  
  // Let's compute grad_Q
  delt = delt0;
  deltQ = Hv(x, delt);
  nablaq = c + deltQ;
  
  // Let's display infos
  pad = '';
  if verbose >0,
    s='   ';
    for i=1:verbose,
      pad = pad + s
    end
    mprintf("\n%s iter    q(delt)  ||nabla q(delt)||  theta\n",pad)
    mprintf("%s %3d  %10.7f    %10.7e   %10.7f\n", pad, niter, (c*delt + 0.5*deltQ*delt), ...
        norm(nablaq), 0)
  end
  
  // Let's compute ||∇q(x)||
  norm2nablaq = nablaq*nablaq';
  pasprecis = norm2nablaq>eps^2;
  
  // Let's compute d_0
  p = -nablaq'
  b = 0
  
  while pasprecis & (niter < MaxIter)
    // d_k
    p = -nablaq' + (b*p);
    pQ = Hv(x, p);
    pQp = pQ*p;
    
    // theta_k
    theta=norm2nablaq/pQp;
    if theta<0.0
      warning("Q not positive ")
    end
    
    // delt_k+1
    delt = delt + theta*p;
    
    // update nabla_q, beta, deltQ
    nablaq = nablaq + theta*pQ;
    nablaq_T = nablaq'
    
    nablaq_Q = Hv(x, nablaq_T)
    b = nablaq_Q*p / pQp
    deltQ = deltQ + theta*pQ;
    
    // update stop conditions    
    norm2nablaq = nablaq*nablaq';
    pasprecis = norm2nablaq>eps^2;
    niter = niter + 1;
    
    if verbose>0,
     mprintf("%s %3d  %10.7f    %10.7e   %10.7f\n", pad, niter, (c*delt + 0.5*deltQ*delt), ...
        norm(nablaq), theta)
    end
    L_f(niter) = c*delt + 0.5*deltQ*delt;
  end
  dN = delt;
  dNQ = deltQ;
endfunction
