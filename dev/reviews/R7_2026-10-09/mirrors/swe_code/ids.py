import sys, os, re
sys.path.insert(0, '/home/claude/lsqvce/dev/docs_src')
os.chdir('/home/claude/lsqvce')
import audit_tree as a
raised = a.raised_ids()
def helpblock(f):
    out=[]; started=False
    for l in open(f, encoding='utf-8'):
        s=l.strip()
        if s.startswith('%'): started=True; out.append(s)
        elif started and s=='' : out.append('')
        elif started: break
    return '\n'.join(out)
def anycomment(f):
    return '\n'.join(l for l in open(f,encoding='utf-8') if l.lstrip().startswith('%'))
files={a._unit_of(f):f for f in a.pglob(a.PKG_GLOB, recursive=True)}
print("NOT in any comment of raising unit:")
for i,us in sorted(raised.items()):
    if not any(i in anycomment(files[u]) for u in us if u in files): print('  ',i,sorted(us))
print("In a comment but NOT in the leading help block:")
for i,us in sorted(raised.items()):
    if any(i in anycomment(files[u]) for u in us if u in files) and not any(i in helpblock(files[u]) for u in us if u in files): print('  ',i,sorted(us))
