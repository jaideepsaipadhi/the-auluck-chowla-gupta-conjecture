/* Rigorous finite check of ACG unimodality for all n <= N.
   Interval arithmetic in x87 long double (64-bit mantissa, exponent range to 2^16383,
   so p(n) ~ e^{1150} for n=2e5 needs NO scaling). lo[] computed with FE_DOWNWARD,
   hi[] with FE_UPWARD. All quantities are sums of nonnegative numbers starting from
   exact v[0]=1, so by induction lo <= exact <= hi.
   Column recurrence: after step k, v[m] = p_{<=k}(m) for m <= N-k; p(n,k) = v[n-k].
   For each n, pairs (k,k+1), 1 <= k <= ceil(n/2), are classified +,-,0 or '?'.
   Output: per-n line "n mode viol" and list of undetermined pairs (n k) to resolve. */
#include <stdio.h>
#include <stdlib.h>
#include <fenv.h>
#include <omp.h>
#pragma STDC FENV_ACCESS ON
typedef long double LD;
int main(int argc,char**argv){
  int N=atoi(argv[1]); FILE*und=fopen(argv[2],"w"); FILE*log=fopen(argv[3],"w");
  int K=(N+1)/2+1;                  /* max k needed */
  LD *lo=calloc(N+1,sizeof(LD)),*hi=calloc(N+1,sizeof(LD));
  LD *plo=calloc(N+1,sizeof(LD)),*phi=calloc(N+1,sizeof(LD));
  LD *blo=calloc(N+1,sizeof(LD));   /* certified lower bound of current max (for mode) */
  int *mode=calloc(N+1,sizeof(int)); char *dec=calloc(N+1,1),*viol=calloc(N+1,1);
  int *nund=calloc(N+1,sizeof(int));
  long long cntp=0,cntm=0,cnt0=0,cntq=0;
  lo[0]=hi[0]=1.0L;
  #pragma omp parallel num_threads(2) reduction(+:cntp,cntm,cnt0,cntq)
  { int t=omp_get_thread_num();
    fesetround(t==0?FE_DOWNWARD:FE_UPWARD);
    volatile LD *v=(t==0)?lo:hi;
    for(int k=1;k<=K;k++){
      int M=N-k;                       /* v valid up to N-k after step k */
      for(int m=k;m<=M;m++) v[m]+=v[m-k];
      #pragma omp barrier
      int n0=k+(N-k+1)*t/2, n1=k+(N-k+1)*(t+1)/2;  /* split n in [k,N] */
      for(int n=n0;n<n1;n++){
        int c=(n+1)/2;                 /* ceil(n/2) */
        if(k>c+1) continue;
        LD al=lo[n-k],ah=hi[n-k];
        if(k>=2){ /* pair (k-1,k), k-1<=c */
          char s;
          if(al>phi[n]) s='+'; else if(ah<plo[n]) s='-';
          else if(al==ah && plo[n]==phi[n] && al==plo[n]) s='0';
          else s='?';
          if(s=='+'){cntp++; if(dec[n]) viol[n]=1;}
          else if(s=='-'){cntm++; dec[n]=1;}
          else if(s=='0') cnt0++;
          else {cntq++; nund[n]++;
            #pragma omp critical
            fprintf(und,"%d %d\n",n,k-1);}
          if(s=='+'){ if(al>blo[n]){blo[n]=al;mode[n]=k;} }
        } else { blo[n]=al; mode[n]=1; }
        plo[n]=al; phi[n]=ah;
      }
      #pragma omp barrier
    }
  }
  int nviol=0;
  for(int n=1;n<=N;n++){ if(viol[n])nviol++;
    fprintf(log,"%d %d %d %d\n",n,mode[n],viol[n],nund[n]); }
  fprintf(stderr,"N=%d pairs: +%lld -%lld 0%lld ?%lld ; float-certified violations: %d\n",N,cntp,cntm,cnt0,cntq,nviol);
  fclose(und); fclose(log); return 0;
}
