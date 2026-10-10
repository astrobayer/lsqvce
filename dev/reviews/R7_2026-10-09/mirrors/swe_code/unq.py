import os,re
root='/home/claude/lsqvce/+lsqvce'
pk={}
for d,_,fs in os.walk(root):
    for f in fs:
        if f.endswith('.m') and f!='Contents.m':
            pk.setdefault(os.path.basename(d),set()).add(f[:-2])
allunits={u for s in pk.values() for u in s}
priv=pk.get('private',set())
for d,_,fs in os.walk(root):
    for f in fs:
        if not f.endswith('.m') or f=='Contents.m': continue
        p=os.path.join(d,f); src=open(p,encoding='utf-8',errors='replace').read()
        locs=set(re.findall(r'^\s*function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?(\w+)',src,re.M))
        privhere=set(os.listdir(os.path.join(d,'private'))) if os.path.isdir(os.path.join(d,'private')) else set()
        privhere={x[:-2] for x in privhere}
        if os.path.basename(d)=='private': privhere|={x[:-2] for x in os.listdir(d)}
        for i,l in enumerate(src.split('\n'),1):
            code=l.split('%')[0]
            code=re.sub(r"'[^']*'",'',code); code=re.sub(r'"[^"]*"','',code)
            for m in re.finditer(r'(?<![\.\w@])([A-Za-z]\w*)\s*\(',code):
                n=m.group(1)
                if n in allunits and n not in locs and n not in privhere and n!=f[:-2]:
                    print(os.path.relpath(p,root),i,n,'|',l.strip()[:120])
