Nwind /* number of windows */
Niter /* maximum number of iterations */
nt[i] /* ni */
ebf[k] /* e+βfk */
ebf2[k] /* buffer for e−βfk */
ebw[i][l][k] /* e−βWk(Ri,l) */
fact[k] /* nk × e+βfk */

for(k=1;k<=Nwind;k=k+1)
{
    fact[k]=nt[k]*ebf[k];
}

for(n=1;n<=Niter;n=n+1) /* start the iteration */
{                       
    for(k=1;k<=Nwind;k=k+1) /* loop over windows */
    {                      
        ebfk=zero;
        for(i=1;i<=Nwind;i=i+1) /* loop over windows */
        {                      
            for(l=1;l<=nt[i];l=l+1) /* loop over times series */
            {                          
                bottom=zero;
                for(j=1;j<=Nwind;j=j+1) /* most inner loop */
                {                          
                    bottom=bottom+ebw[i][l][j]*fact[j]; /* compute the denominator */
                }                                    
                ebfk=ebfk+ebw[i][l][k]/bottom; /* sum the integrand */
            }                               
        }                                
        ebf2[k]=ebfk;  /* save in buffer */
        ebf[k]=1/(ebf[1]*ebfk);  /* shift the new free energy */
        fact[k]=nt[k]*ebf[k];   /* replacement */
    }              
    conv=1; /* convergence flag */
    for(k=1;k<=Nwind;k=k+1)
    {                   
        delta=ABS(kbt*log(ebf[k]/ebf2[k])); /* compare new and old free energies */
        if(delta>=tol) conv=0;    /* test the convergence */
        ebf[k]=ebf2[1]/ebf2[k];   /* shift the free energies */
    }
    if(conv>0) break;         /* if convergence, exit the loop */
}                             /* end of the iteration loop */

