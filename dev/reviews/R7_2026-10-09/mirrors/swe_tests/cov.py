# mirror of audit check 16's unit x category derivation, with and without the never-gated test_speed
import re,glob,os,sys
os.chdir('/home/claude/lsqvce')
decl_re=re.compile(r'^\s*% categories:\s*([a-z, ]+?)\s*(?:\|\s*tags:\s*([a-z, ]+?)\s*)?$')
fn_re=re.compile(r'function\s+(?:\w+\s*=\s*)?(test\w+)\s*\(')
def cov(exclude=()):
    uc={}
    for f in sorted(glob.glob('dev/tests/test_*.m')):
        if os.path.basename(f)[:-2] in exclude: continue
        L=open(f,encoding='utf-8',errors='replace').read().split('\n')
        st=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn_re.match(l)) and not re.match(r'function\s+tests\s*=',l)]
        for k,(i,n) in enumerate(st):
            e=st[k+1][0] if k+1<len(st) else len(L)
            m=decl_re.match(L[i+1]); cats={c.strip() for c in m.group(1).split(',') if c.strip()}
            body='\n'.join(L[i+1:e])
            for u in set(re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w+)',body)):
                uc.setdefault(u,set()).update(cats)
    return uc
a=cov(); b=cov(('test_speed',))
for u in sorted(a):
    d=a[u]-b.get(u,set())
    if d: print('%-28s categories reached ONLY via test_speed (never gated): %s'%(u,sorted(d)))
# exemptions vs coverage: stale exemption rows
ex={}
for row in open('dev/tests/category_exemptions.csv',encoding='utf-8'):
    row=row.strip()
    if not row or row.startswith('#') or row.startswith('unit,'): continue
    p=row.split(',',2); ex[(p[0],p[1])]=p[2]
stale=[(u,c) for (u,c) in ex if c in a.get(u,set())]
print('exemption rows:',len(ex),' rows whose category IS covered (stale):',stale)
# category reached only through a unit name appearing in comment/string? count mentions in comment-only lines

# coverage when only CODE (comments stripped) mentions count
sys.path.insert(0,'dev/docs_src')
def strip(l):
    # crude: drop from first % not in a string
    out=[];q=None
    for i,ch in enumerate(l):
        if q:
            out.append(ch)
            if ch==q: q=None
            continue
        if ch=='%': break
        if ch=='"' : q='"'
        elif ch=="'" and not (out and re.match(r"[\w\)\]\}\.']",out[-1])): q="'"
        out.append(ch)
    return ''.join(out)
def cov2():
    uc={}
    for f in sorted(glob.glob('dev/tests/test_*.m')):
        L=open(f,encoding='utf-8',errors='replace').read().split('\n')
        st=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn_re.match(l)) and not re.match(r'function\s+tests\s*=',l)]
        for k,(i,n) in enumerate(st):
            e=st[k+1][0] if k+1<len(st) else len(L)
            m=decl_re.match(L[i+1]); cats={c.strip() for c in m.group(1).split(',') if c.strip()}
            body='\n'.join(strip(x) for x in L[i+1:e])
            for u in set(re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w+)',body)):
                uc.setdefault(u,set()).update(cats)
    return uc
c=cov2()
tot=0
for u in sorted(a):
    d=a[u]-c.get(u,set())
    if d: tot+=len(d); print('%-28s credited ONLY by comment mentions: %s'%(u,sorted(d)))
print('cells credited only by comments:',tot)
stale2=[(u,cc) for (u,cc) in ex if cc in c.get(u,set())]
print('stale exemption rows against CODE-only coverage:',len(stale2),stale2)

# honest matrix: code-only mentions AND only gated suites (test_speed excluded)
def cov3():
    uc={}
    for f in sorted(glob.glob('dev/tests/test_*.m')):
        if f.endswith('test_speed.m'): continue
        L=open(f,encoding='utf-8',errors='replace').read().split('\n')
        st=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn_re.match(l)) and not re.match(r'function\s+tests\s*=',l)]
        for k,(i,n) in enumerate(st):
            e=st[k+1][0] if k+1<len(st) else len(L)
            m=decl_re.match(L[i+1]); cats={c.strip() for c in m.group(1).split(',') if c.strip()}
            body='\n'.join(strip(x) for x in L[i+1:e])
            for u in set(re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w+)',body)):
                uc.setdefault(u,set()).update(cats)
    return uc
h=cov3()
units=[]
for f in sorted(glob.glob('+lsqvce/**/*.m',recursive=True)):
    if '/private/' in f: continue
    u=(f.replace('+lsqvce/','').replace('.m','').replace('+util/','util.').replace('+engine/','engine.').replace('+kernel/','kernel.'))
    if not u.endswith('Contents'): units.append(u)
C={'contract','reference','precision','speed','robustness','vectorisation','metamorphic'}
print('\nHONEST MATRIX (code-only, gated suites): gaps WITHOUT an exemption row')
n=0
for u in units:
    for c in sorted(C-h.get(u,set())):
        if (u,c) not in ex: n+=1; print('  %-28s %s'%(u,c))
print('total un-exempted gaps:',n)
