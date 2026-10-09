import re,glob,os
os.chdir('/home/claude/lsqvce/dev/tests')
fn=re.compile(r'^function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?(\w+)\s*\(')
for f in sorted(glob.glob('test_*.m')):
    L=open(f,encoding='utf-8').read().split('\n')
    st=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn.match(l))]
    for k,(i,n) in enumerate(st):
        if not n.lower().startswith('test') or n==os.path.basename(f)[:-2]: continue
        e=st[k+1][0] if k+1<len(st) else len(L)
        body='\n'.join(x.split('%')[0] for x in L[i+1:e])
        if not re.search(r'verify\w*\(|assert\w*\(|fatalAssert|assume\w*\(',body):
            print('%s:%d %s has no verify/assert call'%(f,i+1,n))
