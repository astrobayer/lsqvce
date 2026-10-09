import os,re,sys
root='/home/claude/lsqvce/+lsqvce'
units=set()
for d,_,fs in os.walk(root):
    for f in fs:
        if f.endswith('.m'):
            rel=os.path.relpath(os.path.join(d,f),root)
            parts=[p.lstrip('+') for p in rel[:-2].split(os.sep) if p!='private']
            units.add('lsqvce.'+'.'.join(parts))
# methods of classes are also fine: allow lsqvce.X.Class.method
pat=re.compile(r'\blsqvce\.(?:(engine|kernel|util)\.)?([A-Za-z]\w*)',re.I)
bad={}
for d,_,fs in os.walk(root):
    for f in fs:
        if not f.endswith('.m'): continue
        p=os.path.join(d,f)
        for i,l in enumerate(open(p,encoding='utf-8',errors='replace'),1):
            for m in pat.finditer(l):
                tier=(m.group(1) or '').lower(); name=m.group(2)
                cand=('lsqvce.'+tier+'.'+name) if tier else 'lsqvce.'+name
                ok=any(u.lower()==cand.lower() for u in units)
                if not ok and not tier and name.lower() in ('engine','kernel','util'): ok=True
                if not ok:
                    bad.setdefault(cand.lower(),[]).append('%s:%d'%(os.path.relpath(p,root),i))
for k,v in sorted(bad.items()):
    print(k,len(v),v[:6])
