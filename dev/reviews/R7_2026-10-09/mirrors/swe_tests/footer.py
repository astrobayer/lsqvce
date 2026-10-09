import os,re,importlib.util,io,contextlib
os.chdir('/home/claude/lsqvce')
spec=importlib.util.spec_from_file_location('a','dev/docs_src/audit_tree.py'); a=importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
raised=a.raised_ids()
import glob
helpblk={}; allc={}
for f in sorted(glob.glob('+lsqvce/**/*.m',recursive=True)):
    u=a._unit_of(f)
    L=open(f,encoding='utf-8').read().split('\n')
    # help block = first contiguous comment block after the function/classdef line, up to 'See also'
    txt=[];started=False
    for l in L:
        s=l.strip()
        if not started:
            if s.startswith('%'): started=True
            else: continue
        if not s.startswith('%'): break
        if re.match(r'%\s*See also',s): break
        txt.append(s)
    helpblk[u]='\n'.join(txt)
    allc[u]='\n'.join(l for l in L if l.lstrip().startswith('%'))
aud=0; strict=0; foot=[]
for i,us in sorted(raised.items()):
    inany=any(i in allc.get(u,'') for u in us)
    inhelp=any(i in helpblk.get(u,'') for u in us)
    if not inany: aud+=1
    if not inhelp: strict+=1
    if inany and not inhelp: foot.append((i,','.join(sorted(us))))
print('undocumented per check 20 (any comment line):',aud)
print('undocumented if only the help block before "See also" counts:',strict)
print('ids "documented" only in the footer / below See also / in code comments:',len(foot))
for x in foot: print('  ',x)
