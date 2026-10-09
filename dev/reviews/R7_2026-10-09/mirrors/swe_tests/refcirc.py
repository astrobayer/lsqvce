# 'reference'-tagged points: does the body contain any non-lsqvce oracle (literal numbers, dense algebra, closed forms)?
import re,glob,os
os.chdir('/home/claude/lsqvce/dev/tests')
fn=re.compile(r'^function\s+(?:\w+\s*=\s*)?(test\w+)\s*\(')
for f in sorted(glob.glob('test_*.m')):
    L=open(f,encoding='utf-8').read().split('\n')
    st=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn.match(l))]
    for k,(i,n) in enumerate(st):
        e=st[k+1][0] if k+1<len(st) else len(L)
        if 'reference' not in L[i+1]: continue
        body='\n'.join(x.split('%')[0] for x in L[i+1:e])
        oracle=re.search(r'\\|inv\(|chol\(|toeplitz\(|eig\(|qr\(|\d\.\d{3,}|filter\(|pinv\(|\bnorminv|\bchi2|fft\(|legendre|cos\(|sin\(|exact|closed|Ref\b|truth|xTrue|planted|O\d',body)
        if not oracle:
            print('%s:%d %s  <- no visible external oracle in body'%(f,i+1,n))
