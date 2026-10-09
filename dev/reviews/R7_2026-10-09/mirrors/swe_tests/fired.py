import csv,glob,re,os
os.chdir('/home/claude/lsqvce')
rows=[r for r in csv.reader(l for l in open('dev/docs_src/error_ids.csv') if not l.startswith('#'))][1:]
code=''
for f in glob.glob('dev/tests/**/*.m',recursive=True):
    code+='\n'.join(re.sub(r'%.*$','',l) for l in open(f,encoding='utf-8').read().split('\n'))+'\n'
anym=''.join(open(f,encoding='utf-8').read() for f in glob.glob('dev/tests/**/*.m',recursive=True))
inany=sum(1 for r in rows if r[0] in anym); incode=sum(1 for r in rows if r[0] in code)
print('catalogued',len(rows),'| id string anywhere in dev/tests:',inany,'| in test CODE (comments stripped):',incode)
print('comment-only mentions:',[r[0] for r in rows if r[0] in anym and r[0] not in code])
for k in ('error','warning','mixed'):
    sub=[r for r in rows if r[2]==k]; print(k,len(sub),'unfired',sum(1 for r in sub if r[0] not in code))
