import os,re
root='/home/claude/lsqvce/+lsqvce'
rows=[]
for d,_,fs in os.walk(root):
    for f in fs:
        if not f.endswith('.m'): continue
        L=open(os.path.join(d,f),encoding='utf-8',errors='replace').read().split('\n')
        starts=[i for i,l in enumerate(L) if re.match(r'^\s*function\b',l)]
        starts.append(len(L))
        for a,b in zip(starts,starts[1:]):
            m=re.search(r'function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?(\w+)',L[a])
            code=sum(1 for l in L[a:b] if l.strip() and not l.strip().startswith('%'))
            rows.append((b-a,code,os.path.relpath(os.path.join(d,f),root),m.group(1) if m else '?'))
rows.sort(reverse=True)
for r in rows[:14]: print(r)
long=sum(1 for l in open('/dev/null'))
