# tolerances whose tol: comment claims eps/cond*eps-level provenance but whose bound is >= 1e-6 relative/absolute
import re,glob,os,collections
os.chdir('/home/claude/lsqvce/dev/tests')
tolre=re.compile(r"'(Abs|Rel)Tol'\s*,\s*([^,\)]+)")
epsclaim=re.compile(r'(cond\s*\*\s*eps|eps-level|~\s*N\*eps|machine|eps\*cond|roundoff|rounding)',re.I)
cnt=collections.Counter(); hits=[]
for f in sorted(glob.glob('test_*.m')):
    if f=='test_speed.m': continue
    L=open(f,encoding='utf-8').read().split('\n')
    for i,l in enumerate(L):
        for m in tolre.finditer(l):
            try: v=eval(m.group(2).replace('^','**').replace('eps','2.22e-16'),{})
            except Exception: continue
            ctx='\n'.join(L[max(0,i-3):i+1])
            tl=[x for x in L[max(0,i-3):i+1] if 'tol:' in x]
            for t in tl: cnt[t.strip()]+=1
            if isinstance(v,(int,float)) and v>=1e-6 and tl and epsclaim.search(' '.join(tl)) and not re.search(r'MEASURED|measured|Monte|scatter|sigma|probe',' '.join(tl)):
                hits.append((f,i+1,m.group(1),v,tl[-1].strip()[:140]))
for h in hits: print('%s:%d %sTol=%g  <- %s'%h)
print('\nmost-copied tol: comments')
for t,c in cnt.most_common(6): print(c,t[:150])
print('distinct tol comments:',len(cnt),' total sites:',sum(cnt.values()))
