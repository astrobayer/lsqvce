exec(open('/tmp/claude-0/-home-claude/de66562f-20a2-5cf1-9c5b-d39655ca4eeb/scratchpad/mirrors/swe_tests/inject.py').read().split("T='dev/tests")[0])
def sub(p,a,b):
    s=open(p,encoding='utf-8').read(); assert a in s, (p,a); open(p,'w',encoding='utf-8').write(s.replace(a,b,1))
P='+lsqvce/+util/getfielddef.m'; T='dev/tests/test_signalCore.m'
case('15 wrong H1 line in a PACKAGE file (control)',A.check_h1_names,
     lambda d: sub(os.path.join(d,P),'%LSQVCE.UTIL.GETFIELDDEF','%LSQVCE.UTIL.WRONGNAME'))
case('15 wrong H1 line in a TEST file',A.check_h1_names,
     lambda d: sub(os.path.join(d,T),'%TEST_SIGNALCORE','%TEST_WRONGNAME'))
