import re,glob
dq=0; lt=0; lt_no=0; ex=[]
for f in sorted(glob.glob('/home/claude/lsqvce/dev/tests/test_*.m')):
    L=open(f).read().split('\n')
    for i,l in enumerate(L):
        code=l.split('%')[0]
        if re.search(r'"(Abs|Rel)Tol"',code): dq+=1
        m=re.search(r'verify(LessThan|LessThanOrEqual|GreaterThan|GreaterThanOrEqual)\(\s*(tc\s*,\s*)?[^,]*(abs|norm|max)\([^;]*,\s*([0-9.]+e-?[0-9]+|0\.[0-9]+)',code)
        if m:
            lt+=1
            ctx='\n'.join(L[max(0,i-3):i+1])
            if 'tol:' not in ctx:
                lt_no+=1; ex.append('%s:%d: %s'%(f.split('/')[-1],i+1,l.strip()[:110]))
print('double-quoted AbsTol/RelTol sites (invisible to check 16d):',dq)
print('verifyLess/Greater bounds on abs/norm/max with a literal tolerance:',lt,' without any tol: comment within 3 lines:',lt_no)
for e in ex[:12]: print('  ',e)
