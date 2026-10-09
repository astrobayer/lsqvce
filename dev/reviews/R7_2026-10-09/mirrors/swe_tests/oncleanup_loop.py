# Mirror of test_reliability.m:371-379: inside a loop, `w = warning('off','all'); c = onCleanup(@() warning(w));`
# Documented MATLAB semantics assumed: warning('off','all') RETURNS the state before the change;
# assigning a new onCleanup to c destroys the previous object AFTER the RHS is built, running its restore.
def run(branches):
    state='O'; c=None; during=[]
    for b in branches:
        w=state; state='ALLOFF'               # w = warning('off','all')
        new=w                                  # c_k captures w_k
        if c is not None: state=c              # old c destroyed -> warning(w_{k-1})
        c=new
        during.append((b,state))               # state while addSystem/solve run
    state=c                                    # function exit: last c destroyed
    return during,state
for br in (['none','tikhonov','tsvd','lm','auto'], ['none','tikhonov','tsvd','lm','auto','gcv']):
    d,final=run(br); print(len(br),'branches; state during each body:',d,' state after the test:',final)
