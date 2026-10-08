/* ws_fc_exact.c -- floating-point-free verifier for ACG unimodality, n <= N.
   Exact multi-limb unsigned integers (64-bit limbs, per-entry length).
   Column recurrence: after step k, v[m] = p_{<=k}(m) = p(m+k, k).
   At step k, for m ascending: v[m] += v[m-k] (m >= k). Immediately after the
   update of v[m], v[m+1] still holds p_{<=k-1}(m+1) = p(n, k-1) with n = m+k,
   so p(n,k) vs p(n,k-1) is compared exactly with no extra storage.
   Pairs (k-1,k) are checked for every n <= N with k-1 <= ceil(n/2), k >= 2.
   Output log: "n mode viol nund" (nund always 0), same format as ws_fc_check.
   No floating point anywhere (compile check: grep for float/double finds none). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64;
static int L; static u64 *V; static unsigned short *len;
#define E(m) (V+(size_t)(m)*L)
static inline void add(int d,int s){ /* v[d]+=v[s] */
  u64 *a=E(d),*b=E(s); int la=len[d],lb=len[s],n=la>lb?la:lb; unsigned char c=0;
  for(int i=0;i<n;i++){ u64 x=i<la?a[i]:0, y=i<lb?b[i]:0, r;
    unsigned char c1=__builtin_add_overflow(x,y,&r); unsigned char c2=__builtin_add_overflow(r,(u64)c,&r);
    a[i]=r; c=c1|c2; }
  if(c){ if(n>=L){fprintf(stderr,"limb overflow\n");exit(2);} a[n]=1; n++; }
  len[d]=n; }
static inline int cmp(int x,int y){
  if(len[x]!=len[y]) return len[x]>len[y]?1:-1;
  u64 *a=E(x),*b=E(y); for(int i=len[x]-1;i>=0;i--) if(a[i]!=b[i]) return a[i]>b[i]?1:-1;
  return 0; }
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: %s N log\n",argv[0]);return 1;}
  int N=atoi(argv[1]); FILE*log=fopen(argv[2],"w");
  /* limb budget: p(N) < exp(pi*sqrt(2N/3)) <= 2^(1.4427*pi*sqrt(2N/3)); +2 slack, overflow is checked */
  { long b=1; long s=1; while(s*s< (long)(2L*N/3+1)) s++; b=(long)(4533*s/1000)+2; L=(int)(b/64+2); }
  V=calloc((size_t)(N+2)*L,8); len=calloc(N+2,sizeof *len);
  int *mode=calloc(N+1,sizeof(int)); char *dec=calloc(N+1,1),*viol=calloc(N+1,1);
  long long cp=0,cm=0,c0=0;
  V[0]=1; len[0]=1; for(int n=1;n<=N;n++) mode[n]=1;   /* k=1: v[m]=1 for all m */
  for(int m=1;m<=N;m++){V[(size_t)m*L]=1;len[m]=1;}
  int K=(N+1)/2+1;
#define PAIR(m) do{ int n=(m)+k; if(n<=N && 2*(k-1)<=n+1){ int s=cmp((m),(m)+1); \
      if(s>0){cp++; if(dec[n]) viol[n]=1; else mode[n]=k;} else if(s<0){cm++;dec[n]=1;} else c0++; } }while(0)
  for(int k=2;k<=K;k++){
    for(int m=(k>=3?k-3:0); m<k && m<=N-k; m++) PAIR(m);   /* m<k: v[m] unchanged at step k */
    for(int m=k;m<=N-k;m++){ add(m,m-k); PAIR(m); }
  }
  int nv=0; for(int n=1;n<=N;n++){ nv+=viol[n]; fprintf(log,"%d %d %d 0\n",n,mode[n],viol[n]); }
  fclose(log);
  fprintf(stderr,"N=%d limbs/slot=%d pairs: +%lld -%lld 0%lld ; exact violations: %d -> %s\n",N,L,cp,cm,c0,nv,nv?"FAIL":"PASS");
  return nv?3:0; }
