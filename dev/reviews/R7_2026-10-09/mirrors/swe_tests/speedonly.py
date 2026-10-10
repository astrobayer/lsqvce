import re,glob,os
os.chdir('/home/claude/lsqvce')
decl_re = re.compile(r'^\s*% categories:\s*([a-z, ]+?)\s*(?:\|\s*tags:\s*([a-z, ]+?)\s*)?$')
fn_re = re.compile(r'function\s+(?:\w+\s*=\s*)?(test\w+)\s*\(')
def matrix(skip=()):
    uc={}
    for f in sorted(glob.glob('dev/tests/test_*.m')):
        if os.path.basename(f) in skip: continue
        L=open(f).read().split('\n')
        st=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn_re.match(l)) and not re.match(r'function\s+tests\s*=',l)]
        for k,(i,n) in enumerate(st):
            end=st[k+1][0] if k+1<len(st) else len(L)
            m=decl_re.match(L[i+1]); cats={c.strip() for c in m.group(1).split(',') if c.strip()}
            for u in set(re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w+)','\n'.join(L[i+1:end]))):
                uc.setdefault(u,set()).update(cats)
    return uc
a=matrix(); b=matrix(skip=('test_speed.m',))
for u in sorted(a):
    d=a[u]-b.get(u,set())
    if d: print(u,'credited ONLY by the never-gated test_speed:',sorted(d))
# also: units credited only via comments? count matches of lsqvce.X in comment lines only
print('---- comment-only credit')
def matrix2():
    uc={}
    for f in sorted(glob.glob('dev/tests/test_*.m')):
        L=open(f).read().split('\n')
        st=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn_re.match(l)) and not re.match(r'function\s+tests\s*=',l)]
        for k,(i,n) in enumerate(st):
            end=st[k+1][0] if k+1<len(st) else len(L)
            m=decl_re.match(L[i+1]); cats={c.strip() for c in m.group(1).split(',') if c.strip()}
            code='\n'.join(l.split('%')[0] for l in L[i+1:end])   # crude: drop comments (and anything after % in strings)
            for u in set(re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w+)',code)):
                uc.setdefault(u,set()).update(cats)
    return uc
c=matrix2()
for u in sorted(a):
    d=a[u]-c.get(u,set())
    if d: print(u,'credited only via a COMMENT/string mention:',sorted(d))
print('---- strict: code lines only, body ends at the NEXT function of any kind, test_speed excluded')
anyfn=re.compile(r'^\s*function\b')
def matrix3():
    uc={}
    for f in sorted(glob.glob('dev/tests/test_*.m')):
        if f.endswith('test_speed.m'): continue
        L=open(f).read().split('\n')
        st=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn_re.match(l)) and not re.match(r'function\s+tests\s*=',l)]
        for (i,n) in st:
            j=i+1
            while j<len(L) and not anyfn.match(L[j]): j+=1
            m=decl_re.match(L[i+1]); cats={c.strip() for c in m.group(1).split(',') if c.strip()}
            code='\n'.join(re.sub(r"%.*$","",l) if not l.lstrip().startswith('%') else '' for l in L[i+1:j])
            for u in set(re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w+)',code)):
                uc.setdefault(u,set()).update(cats)
    return uc
s=matrix3()
units=[]
for f in sorted(glob.glob('+lsqvce/**/*.m',recursive=True)):
    if '/private/' in f: continue
    u=f.replace('+lsqvce/','').replace('.m','').replace('+util/','util.').replace('+engine/','engine.').replace('+kernel/','kernel.')
    if not u.endswith('Contents'): units.append(u)
exempt=set()
for row in open('dev/tests/category_exemptions.csv'):
    row=row.strip()
    if not row or row[0]=='#' or row.startswith('unit,'): continue
    p=row.split(',',2); exempt.add((p[0],p[1]))
tot=0
for u in units:
    lost=sorted(a.get(u,set())-s.get(u,set()))
    unex=[c for c in lost if (u,c) not in exempt]
    if unex: print(' ',u,'loses',unex,'(no exemption row)'); tot+=len(unex)
print('gaps the audit would report under strict parsing:',tot)
print('---- attribution: which of the strict losses are comment-only (body as audit defines it, test_speed included)')
for u in units:
    lost=sorted(a.get(u,set())-c.get(u,set()))
    unex=[x for x in lost if (u,x) not in exempt]
    if unex: print(' COMMENT-ONLY',u,unex)
    lost2=sorted(a.get(u,set())-b.get(u,set()))
    unex2=[x for x in lost2 if (u,x) not in exempt]
    if unex2: print(' TEST_SPEED-ONLY',u,unex2)
print('---- helper-function attribution examples')
for f in sorted(glob.glob('dev/tests/test_*.m')):
    if f.endswith('test_speed.m'): continue
    L=open(f).read().split('\n')
    st=[(i,m.group(1)) for i,l in enumerate(L) if (m:=fn_re.match(l)) and not re.match(r'function\s+tests\s*=',l)]
    for k,(i,n) in enumerate(st):
        end=st[k+1][0] if k+1<len(st) else len(L)
        j=i+1
        while j<len(L) and not anyfn.match(L[j]): j+=1
        cats=decl_re.match(L[i+1]).group(1)
        own=set(re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w+)','\n'.join(L[i+1:j])))
        tail=set(re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w+)','\n'.join(L[j:end])))
        extra=tail-own
        for u in ['engine.BlockAccumulator','engine.synthLoader','engine.adjustBlock']:
            if u in extra and any(c in cats for c in ['vectorisation','speed']):
                hf=re.match(r'\s*function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?(\w+)',L[j]).group(1)
                print(' ',os.path.basename(f),n,'['+cats+'] credits',u,'via helper',hf,'at line',j+1)
