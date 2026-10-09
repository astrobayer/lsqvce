# Mirror of lsqvce.Options: constructor assigns N=V pairs in order; set.Parallel = mergeGroup(DEFAULTS, s)
DEF={'Mode':'serial','DistributeDays':False,'Threshold':8}
def merge(d,s):
    g=dict(d)
    for k,v in s.items():
        assert k in d; g[k]=v
    return g
obj={}
for name,val in [('Verbose',False),('Parallel',{'Mode':'parfor'}),('Parallel',{'DistributeDays':True,'Threshold':2})]:
    obj[name]=merge(DEF,val) if name=='Parallel' else val
P=obj['Parallel']; print('resolved Parallel:',P)
distrib = P['DistributeDays'] and P['Mode'].lower()!='serial' and 3>=P['Threshold']
print('BlockAccumulator distrib (3 days):',distrib)
