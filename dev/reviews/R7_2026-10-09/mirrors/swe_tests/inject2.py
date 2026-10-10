exec(open('/tmp/claude-0/-home-claude/de66562f-20a2-5cf1-9c5b-d39655ca4eeb/scratchpad/mirrors/swe_tests/inject.py').read().split("T='dev/tests")[0])
P='+lsqvce/+util/getfielddef.m'
case('15 wrong H1 line in a PACKAGE file (control, fixed)',A.check_h1_names,
     lambda d: open(os.path.join(d,P),'w').write(open(os.path.join(d,P)).read().replace('%LSQVCE.UTIL.GETFIELDDEF','%LSQVCE.UTIL.WRONGNAME',1)))
case('7  unbalanced end in a PACKAGE file (control)',A.check_balance,
     lambda d: app(os.path.join(d,P),'\nfunction zzUnbalanced()\n if true\nend\n'))
case('2  removed name in PACKAGE code (control)',A.check_removed_names,
     lambda d: app(os.path.join(d,P),'\nfunction zz()\n  x = coalesceOpt(1);\nend\n'))
