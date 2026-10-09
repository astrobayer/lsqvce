import re,glob,os,csv,collections
os.chdir('/home/claude/lsqvce')
rows=[r for r in open('dev/docs_src/error_ids.csv',encoding='utf-8') if r.strip() and not r.startswith('#') and not r.startswith('id,')]
ids=[]
for r in csv.reader(rows):
    ids.append((r[0],r[1],r[2]))
tiers={}
for r in csv.DictReader(open('dev/docs_src/tiers.csv',encoding='utf-8')):
    tiers[r['name']]=r['tier']
tests=''
for f in glob.glob('dev/tests/*.m')+glob.glob('dev/tests/private/*.m'):
    tests+=open(f,encoding='utf-8',errors='replace').read()
# code-only text of tests (strip comments)
def strip_comments(t):
    out=[]
    for l in t.split('\n'):
        s=l.lstrip()
        if s.startswith('%'): continue
        out.append(l)
    return '\n'.join(out)
tc=strip_comments(tests)
unc=collections.defaultdict(list)
for i,u,k in ids:
    if i not in tc:
        units=u.split(';')
        tier=min((tiers.get(x,'?') for x in units), key=lambda t:{'USER':0,'ENGINE':1,'KERNEL':2}.get(t,3))
        unc[tier].append((i,u,k))
n=sum(len(v) for v in unc.values())
print('catalogued ids:',len(ids),' not appearing in any test code line:',n)
for t in ['USER','ENGINE','KERNEL','?']:
    for i,u,k in sorted(unc.get(t,[])): print('%-7s %-8s %-55s %s'%(t,k,i,u))
