import os,re,glob,collections
os.chdir('/home/claude/lsqvce')
units={}
for f in sorted(glob.glob('+lsqvce/**/*.m',recursive=True)):
    if '/private/' in f or f.endswith('Contents.m'): continue
    u=f.replace('+lsqvce/','').replace('.m','').replace('+util/','util.').replace('+engine/','engine.').replace('+kernel/','kernel.')
    units[u]=f
code={}
for f in glob.glob('+lsqvce/**/*.m',recursive=True):
    code[f]='\n'.join(re.sub(r'%.*$','',l) for l in open(f,encoding='utf-8').read().split('\n'))
rows=[l.strip().split(',',2) for l in open('dev/tests/category_exemptions.csv') if l.strip() and not l.startswith('#') and not l.startswith('unit,')]
callerclaim=[r for r in rows if re.search(r'caller|through (its|the) caller',r[2])]
for u in sorted({r[0] for r in callerclaim}):
    short=u.split('.')[-1]
    pat=re.compile(r'\blsqvce\.'+re.escape(u)+r'\b' + (r'|\b'+re.escape(short)+r'\s*\(' if True else ''))
    callers=[f for f,t in code.items() if f!=units.get(u) and re.search(r'\blsqvce\.'+re.escape(u)+r'\b',t)]
    if not callers:
        n=sum(1 for r in callerclaim if r[0]==u)
        print('NO package caller:',u,'-- yet',n,'exemption row(s) argue coverage "through its callers"')
