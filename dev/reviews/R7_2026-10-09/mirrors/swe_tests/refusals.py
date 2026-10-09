import csv,glob,os
os.chdir('/home/claude/lsqvce')
rows=[r for r in csv.reader(l for l in open('dev/docs_src/error_ids.csv') if not l.startswith('#'))][1:]
tests=''.join(open(f,encoding='utf-8').read() for f in glob.glob('dev/tests/**/*.m',recursive=True))
ex=[l.split(',',2) for l in open('dev/tests/category_exemptions.csv') if 'its refusals are named identifiers pinned under contract' in l]
for u,c,r in ex:
    ids=[(i,k) for i,us,k,_ in rows if u in us.split(';') and k!='warning']
    fired=[i for i,k in ids if i in tests]
    print('%-28s errors=%d fired-by-tests=%d  unfired=%s'%(u,len(ids),len(fired),[i.split(':')[-1] for i,k in ids if i not in tests]))
