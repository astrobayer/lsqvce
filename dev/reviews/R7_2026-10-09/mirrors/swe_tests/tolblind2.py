import re,glob,collections
tot=0; no=0; per=collections.Counter()
for f in sorted(glob.glob('/home/claude/lsqvce/dev/tests/test_*.m')):
    L=open(f).read().split('\n')
    for i,l in enumerate(L):
        code=l.split('%')[0]
        if re.search(r'verify(LessThan|LessThanOrEqual|GreaterThan|GreaterThanOrEqual)\b',code) and re.search(r',\s*[0-9]*\.?[0-9]+(e-?[0-9]+)?\s*[,)]',code):
            tot+=1
            if 'tol:' not in '\n'.join(L[max(0,i-3):i+1]): no+=1; per[f.split('/')[-1]]+=1
print('verifyLess/Greater* with a numeric literal bound:',tot,'; with no tol: comment within 3 lines:',no)
print(per.most_common())
