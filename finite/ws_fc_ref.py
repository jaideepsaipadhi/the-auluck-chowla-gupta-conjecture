"""Exact big-integer reference / resolver.
ref N            : exact p(n,k) for n<=N; prints per-n 'n mode viol 0' (mode = first argmax over k<=ceil(n/2)+1)
                   and totals of +,-,0 over pairs (k,k+1), k<=ceil(n/2).  Also checks identity p(n,k)=p(n-k), k>=n/2.
resolve FILE     : exact sign of p(n,k+1)-p(n,k) for each 'n k' line."""
import sys
def table(N):
    # P[n][k] = p(n,k): p(n,k)=p(n-1,k-1)+p(n-k,k)
    P=[[0]*(N+1) for _ in range(N+1)]; P[0][0]=1
    for n in range(1,N+1):
        for k in range(1,n+1): P[n][k]=P[n-1][k-1]+P[n-k][k]
    return P
if sys.argv[1]=="ref":
    N=int(sys.argv[2]); P=table(N); pp=[sum(r) for r in P]; tp=tm=t0=0; out=[]
    for n in range(1,N+1):
        c=(n+1)//2; row=P[n]; dec=False; viol=0
        for k in range(1,c+1):
            d=row[k+1]-row[k]
            if d>0: tp+=1; viol|=dec
            elif d<0: tm+=1; dec=True
            else: t0+=1
        for k in range(1,n+1):
            if 2*k>=n: assert row[k]==pp[n-k]
        sub=row[1:c+2]; mode=sub.index(max(sub))+1
        out.append(f"{n} {mode} {int(viol)} 0")
        assert all(row[k]>=row[k+1] for k in range(c,n))  # full row unimodal tail
    open(sys.argv[3],"w").write("\n".join(out)+"\n")
    print(f"exact N={N}: +{tp} -{tm} 0{t0}")
else:
    pairs=[tuple(map(int,l.split())) for l in open(sys.argv[2]) if l.strip()]
    if not pairs: print("nothing to resolve"); sys.exit()
    N=max(n for n,_ in pairs); P=table(N)
    for n,k in pairs: print(n,k,(P[n][k+1]>P[n][k])-(P[n][k+1]<P[n][k]))
