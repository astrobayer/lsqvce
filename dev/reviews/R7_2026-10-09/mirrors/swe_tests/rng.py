import re,glob,os
os.chdir('/home/claude/lsqvce/dev/tests')
fn=re.compile(r'^function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?(\w+)\s*\(')
draw=re.compile(r'\b(randn|rand|randi|randperm)\s*\(\s*(?!rs\b|rs2\b|stream\b|RandStream\b|s\b|r\b|g\b)')
for f in sorted(glob.glob('test_*.m')):
    L=open(f,encoding='utf-8').read().split('\n')
    funcs=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn.match(l))]
    for k,(i,n) in enumerate(funcs):
        e=funcs[k+1][0] if k+1<len(funcs) else len(L)
        body=[x.split('%')[0] for x in L[i:e]]
        seeded=False
        for j,b in enumerate(body):
            if re.search(r'\brng\s*\(',b): seeded=True
            if draw.search(b) and not seeded:
                print('%s:%d %s: unseeded global draw: %s'%(f,i+j+1,n,b.strip()[:80])); break
