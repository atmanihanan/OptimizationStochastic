
function [dN, dNQ, niter, L_f] = GC_TR(c, Hv,x, eps, MaxIter, Delta,verbose)
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
    [bidon,dim] = size(c);
    n = length(x)
    delt = zeros(n, 1);
    deltQ =Hv(x,delt);
    nablaq = c + deltQ;
    p= -nablaq';
    b = 0;
   
    // Let's display infos
    pad = '';
    if verbose > 0,
        s='   ';
        for i=1:verbose,
            pad = pad + s;
        end
        mprintf("%s iter    q(delt)  ||nabla q(delt)||  theta\n",pad);
        mprintf("%s %3d  %10.7f    %10.7e   %10.7f\n", pad, niter, (c*delt + 0.5*deltQ*delt), ...
            norm(nablaq), 0);
    end
    
    norm2nablaq = nablaq*nablaq';
    pasprecis = norm2nablaq>eps^2;
    sortie = %f
    
    // Boucle
    while pasprecis & (niter < MaxIter) & ~sortie
        // count iteration
        niter = niter + 1;
        
        // refresh variable
        pQ =Hv(x,p);
        pQp = pQ*p;
       
        // compute theta_max
        a = p'*p;
        b=2*delt'*p ;
        c = delt'*delt - Delta^2;
        disc= b^2-4*a*c;
        if disc < 0 then
            sortie = %t
        else
            theta1= (-b-sqrt(disc))/(2*a);
            theta2= (-b+sqrt(disc))/(2*a);
            theta=max(theta1, theta2);
        end
        
        // sortie
        sortie=Hv(x,p)*p <= 0
        
        if ~sortie then
            thetaGC=(-nablaq*p) / pQp;
            
            // BLOC INESISTANT DANS L'ALGO DE BASE
            if thetaGC<0 | thetaGC>theta then
                thetaGC = theta
            end
            // -----------------------------------
            
            sortie=thetaGC>theta;
            if ~sortie then
                delt=delt+thetaGC*p;
                
                // update nablaq & deltQ
                deltQ = deltQ + thetaGC*pQ;
                nablaq = nablaq + thetaGC*pQ;
            
                b = (Hv(x, nablaq')*p) / pQp;
                p= -nablaq' +b*p;
            end
        end
         
        // update stop conditions    
        norm2nablaq = nablaq*nablaq';
        pasprecis = norm2nablaq > eps^2;   
        if verbose > 0,
            mprintf("%s %3d  %10.7f    %10.7e   %10.7f\n", pad, niter, (c*delt + 0.5*deltQ*delt), ...
            norm(nablaq), theta);
        end
        // L_f(niter) = c*delt + 0.5*deltQ*delt;
    end
   
    //disp("sortie")
    /*if sortie then
        disp('sortie happens')
        delt = delt + (thetaGC*p)
    end*/
    dN = delt;
    dNQ = deltQ;
endfunction



