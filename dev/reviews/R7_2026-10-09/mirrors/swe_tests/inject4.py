exec(open('/tmp/claude-0/-home-claude/de66562f-20a2-5cf1-9c5b-d39655ca4eeb/scratchpad/mirrors/swe_tests/inject.py').read().split("T='dev/tests")[0])
P='+lsqvce/+util/getfielddef.m'; T='dev/tests/test_signalCore.m'
def full():
    A.main()
case('FULL AUDIT: help cites id only a TEST raises',full,
     lambda d: (app(os.path.join(d,P),'\n%   lsqvce:zz:ghost  never raised by the package\n'),
                app(os.path.join(d,T),"\nfunction zz(tc)\n error('lsqvce:zz:ghost','x');\nend\n")))
case('FULL AUDIT: duplicate local fn + unbalanced end in a TEST file',full,
     lambda d: app(os.path.join(d,T),'\nfunction makeCfg()\nend\nfunction makeCfg()\nend\nfunction zz(tc)\n if true\nend\n'))
