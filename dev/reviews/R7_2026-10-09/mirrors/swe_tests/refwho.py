import re,glob,os,sys
os.chdir('/home/claude/lsqvce/dev/tests')
fn=re.compile(r'^function\s+(?:\w+\s*=\s*)?(test\w+)\s*\(')
want=sys.argv[1]; cat=sys.argv[2]
for f in sorted(glob.glob('test_*.m')):
    L=open(f,encoding='utf-8').read().split('\n')
    st=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn.match(l))]
    for k,(i,n) in enumerate(st):
        e=st[k+1][0] if k+1<len(st) else len(L)
        if cat not in L[i+1]: continue
        body='\n'.join(L[i+1:e])
        if re.search(r'lsqvce\.'+re.escape(want)+r'\b',body): print(f,i+1,n)
