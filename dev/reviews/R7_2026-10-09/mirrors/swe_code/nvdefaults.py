# Documented Name-Value set and defaults vs the arguments block, every package function.
import re, glob
for f in sorted(glob.glob('/home/claude/lsqvce/+lsqvce/**/*.m', recursive=True)):
    src = open(f, encoding='utf-8').read().split('\n')
    # leading help block
    hb = []
    for l in src[1:]:
        if l.strip().startswith('%'): hb.append(l.strip())
        elif l.strip() == '' and hb: continue
        else: break
    help_ = '\n'.join(hb)
    m = None
    for i, l in enumerate(src):
        if re.match(r'^\s*arguments\s*$', l): m = i; break
    if m is None: continue
    args = {}
    for l in src[m+1:]:
        if re.match(r'^\s*end\s*$', l): break
        mm = re.match(r'^\s*opt\.(\w+)\b[^=%]*?(=\s*([^%]+))?', l)
        if mm: args[mm.group(1)] = (mm.group(3) or '').strip()
    if not args: continue
    rel = f.replace('/home/claude/lsqvce/', '')
    for name, dflt in args.items():
        if not re.search(r"['\"]?%s['\"]?" % name, help_, re.I):
            print('%-42s %-18s accepted, NOT in help' % (rel, name))
            continue
        # documented default near the name
        for hl in hb:
            if re.search(r"'%s'" % name, hl, re.I) or re.search(r"^\%%\s+%s\b" % name, hl):
                md = re.search(r'default\s*([-\w.\'"\[\]]+)', hl, re.I)
                if md and dflt:
                    a = dflt.strip().strip('"\'').lower(); h = md.group(1).strip().strip('"\'.').lower()
                    if a != h and not (a.replace('e-0','e-') == h.replace('e-0','e-')):
                        try:
                            if float(a) == float(h): break
                        except Exception: pass
                        print('%-42s %-18s help default %-10s args default %s' % (rel, name, md.group(1), dflt))
                break
