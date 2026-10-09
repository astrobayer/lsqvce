import re,glob,sys
# crude: split files at top-level 'function' lines (column 0) and classdef method 'function' at 8 spaces
out=[]
for f in glob.glob('/home/claude/lsqvce/+lsqvce/**/*.m', recursive=True):
    L=open(f,encoding='utf-8').read().split('\n')
    starts=[i for i,l in enumerate(L) if re.match(r'^(        )?function\b', l)]
    starts.append(len(L))
    for a,b in zip(starts,starts[1:]):
        # count code lines only (exclude comment/blank)
        code=sum(1 for l in L[a:b] if l.strip() and not l.strip().startswith('%'))
        name=re.sub(r'.*?function\s+(\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?','',L[a]).split('(')[0]
        if (b-a)>350: out.append((code,b-a,f.replace('/home/claude/lsqvce/',''),name.strip()))
for o in sorted(out,reverse=True): print('%5d code / %5d total  %s  %s'%o)
