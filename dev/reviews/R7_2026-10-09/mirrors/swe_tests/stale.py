import re,sys,io,contextlib,os
os.chdir('/home/claude/lsqvce')
sys.path.insert(0,'dev/docs_src')
import importlib.util
spec=importlib.util.spec_from_file_location('audit','dev/docs_src/audit_tree.py')
a=importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
# re-derive unit_cat as check 16 does
decl_re = re.compile(r'^\s*% categories:\s*([a-z, ]+?)\s*(?:\|\s*tags:\s*([a-z, ]+?)\s*)?$')
fn_re = re.compile(r'function\s+(?:\w+\s*=\s*)?(test\w+)\s*\(')
unit_cat={}
import glob
for f in sorted(glob.glob('dev/tests/test_*.m')):
    lines=open(f,encoding='utf-8').read().split('\n')
    starts=[(i,m.group(1)) for i,l in enumerate(lines) if (m:=fn_re.match(l)) and not re.match(r'function\s+tests\s*=',l)]
    for k,(i,name) in enumerate(starts):
        end=starts[k+1][0] if k+1<len(starts) else len(lines)
        m=decl_re.match(lines[i+1])
        cats={c.strip() for c in m.group(1).split(',') if c.strip()}
        body='\n'.join(lines[i+1:end])
        for u in set(re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w+)', body)):
            unit_cat.setdefault(u,set()).update(cats)
stale=[]
for row in open('dev/tests/category_exemptions.csv'):
    row=row.strip()
    if not row or row.startswith('#') or row.startswith('unit,'): continue
    u,c,r=row.split(',',2)
    if c in unit_cat.get(u,set()): stale.append((u,c,r[:90]))
print('STALE exemption rows (category now covered):',len(stale))
for s in stale: print(' ',s)
