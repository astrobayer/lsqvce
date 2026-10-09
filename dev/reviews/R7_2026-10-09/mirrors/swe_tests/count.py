import re,sys,glob,os,collections,json
root='/home/claude/lsqvce/dev/tests'
FUNC=re.compile(r'^\s*function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?(\w+)\s*(\(|$)')
setupnames={'setup','teardown','setuponce','teardownonce'}
tot=0;totM=0;cats=collections.Counter();tags=collections.Counter()
rows=[]
nocat=[]
for f in sorted(glob.glob(root+'/test_*.m')):
    L=open(f,encoding='utf-8',errors='replace').read().split('\n')
    inb=False; names=[]; ml=[]; dp=0
    for i,l in enumerate(L):
        s=l.strip()
        if s=='%{': inb=True; continue
        if s=='%}': inb=False; continue
        if inb: continue
        m=FUNC.match(l)
        if not m: continue
        nm=m.group(1)
        # declaredPoints mirror
        if not re.match(r'^function\s+tests\s*=',s) and re.match(r'^function\s+(?:\w+\s*=\s*)?test\w+\s*\(',s): dp+=1
        if nm.lower() in setupnames: continue
        if nm==os.path.basename(f)[:-2]: continue
        if nm.lower().startswith('test') or nm.lower().endswith('test'):
            names.append(nm)
            # categories line
            nxt=L[i+1].strip() if i+1<len(L) else ''
            mm=re.match(r'%\s*categories:\s*([^|]*)(\|\s*tags:\s*(.*))?',nxt)
            if mm:
                for c in re.split(r'[,\s]+',mm.group(1).strip()):
                    if c: cats[c]+=1
                if mm.group(3):
                    for t in re.split(r'[,\s]+',mm.group(3).strip()):
                        if t: tags[t]+=1
            else: nocat.append((os.path.basename(f),nm,nxt[:60]))
    base=os.path.basename(f)[:-2]
    rows.append((base,len(names),dp))
    tot+=len(names)
for r in rows: print('%-22s matlabRule=%3d declaredPoints=%3d %s'%(r[0],r[1],r[2],'' if r[1]==r[2] else '<-- DIFF'))
print('TOTAL',tot,'declaredPoints total',sum(r[2] for r in rows))
print('excl test_speed:',tot-[r[1] for r in rows if r[0]=='test_speed'][0])
print('cats',dict(cats));print('tags',dict(tags))
print('no categories line:',nocat)
