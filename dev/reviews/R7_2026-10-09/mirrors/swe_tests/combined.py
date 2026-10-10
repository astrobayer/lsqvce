import re,glob,os
os.chdir('/home/claude/lsqvce')
decl_re = re.compile(r'^\s*% categories:\s*([a-z, ]+?)\s*(?:\|\s*tags:\s*([a-z, ]+?)\s*)?$')
fn_re = re.compile(r'function\s+(?:\w+\s*=\s*)?(test\w+)\s*\(')
def strip(l):
    # drop comment text: whole-line comments and trailing % outside quotes (crude quote tracking)
    if l.lstrip().startswith('%'): return ''
    out='';q=None
    for i,ch in enumerate(l):
        if q:
            out+=ch
            if ch==q: q=None
            continue
        if ch=='%': break
        if ch=='"' or (ch=="'" and not (i>0 and (l[i-1].isalnum() or l[i-1] in ")]}'_."))): q=ch
        out+=ch
    return out
def matrix(nocomment, nospeed):
    uc={}
    for f in sorted(glob.glob('dev/tests/test_*.m')):
        if nospeed and f.endswith('test_speed.m'): continue
        L=open(f).read().split('\n')
        st=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn_re.match(l)) and not re.match(r'function\s+tests\s*=',l)]
        for k,(i,n) in enumerate(st):
            end=st[k+1][0] if k+1<len(st) else len(L)
            cats={c.strip() for c in decl_re.match(L[i+1]).group(1).split(',') if c.strip()}
            seg=L[i+1:end]
            body='\n'.join(strip(l) for l in seg) if nocomment else '\n'.join(seg)
            for u in set(re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w+)',body)):
                uc.setdefault(u,set()).update(cats)
    return uc
units=[]
for f in sorted(glob.glob('+lsqvce/**/*.m',recursive=True)):
    if '/private/' in f: continue
    u=f.replace('+lsqvce/','').replace('.m','').replace('+util/','util.').replace('+engine/','engine.').replace('+kernel/','kernel.')
    if not u.endswith('Contents'): units.append(u)
ex=set()
for row in open('dev/tests/category_exemptions.csv'):
    row=row.strip()
    if not row or row[0]=='#' or row.startswith('unit,'): continue
    p=row.split(',',2); ex.add((p[0],p[1]))
cats={'contract','reference','precision','speed','robustness','vectorisation','metamorphic'}
def gaps(m): return sorted((u,c) for u in units for c in cats if c not in m.get(u,set()) and (u,c) not in ex)
print('audit as written:',len(gaps(matrix(False,False))))
g=gaps(matrix(True,True)); print('comments stripped + never-gated test_speed excluded:',len(g))
for x in g: print('  ',x)
