import sys, os, re
sys.path.insert(0, '/home/claude/lsqvce/dev/docs_src'); os.chdir('/home/claude/lsqvce')
import audit_tree as a
raised = a.raised_ids()
allr = set(raised)
for f in sorted(a.pglob(a.PKG_GLOB, recursive=True)):
    txt = open(f, encoding='utf-8').read()
    hb = '\n'.join(l for l in txt.split('\n') if l.lstrip().startswith('%'))
    hb = re.split(r'See also', hb)[0]
    for i in sorted(set(re.findall(r'lsqvce:[A-Za-z]+:[A-Za-z]+(?::[A-Za-z]+)?', hb))):
        if i not in allr:
            print('%-45s documents %s -- raised NOWHERE' % (f, i))
