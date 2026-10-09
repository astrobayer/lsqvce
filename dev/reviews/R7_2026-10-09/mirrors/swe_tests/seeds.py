import re,glob
FN=re.compile(r'^\s*function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?(\w+)\s*\(')
for f in sorted(glob.glob('/home/claude/lsqvce/dev/tests/test_*.m')+glob.glob('/home/claude/lsqvce/dev/tests/fixture_*.m')):
    L=open(f).read().split('\n')
    starts=[(i,m.group(1)) for i,l in enumerate(L) if (m:=FN.match(l))]
    for k,(i,n) in enumerate(starts):
        end=starts[k+1][0] if k+1<len(starts) else len(L)
        body=[l.split('%')[0] for l in L[i:end]]
        txt='\n'.join(body)
        glob_rand=re.findall(r'\b(randn|rand|randi|randperm)\(\s*(?!rs\b|rs2\b|s\b|stream|str\b|r\b|RS)([^)]*)',txt)
        glob_rand=[g for g in glob_rand if not re.match(r'\s*(rs\w*|s\w*|stream\w*)\s*,',g[1])]
        seeded=re.search(r'\brng\(',txt)
        if glob_rand and not seeded:
            print(f.split('/')[-1],n,'global-RNG draws without rng() in this function:',[g[0]+'('+g[1][:20] for g in glob_rand][:3])
        if seeded and 'rng(' in txt and not re.search(r'=\s*rng\b|rng\(\s*(s0|prev|old|before)',txt) and n.lower().startswith('test'):
            pass
