# Fault-inject audit checks on a SCRATCH COPY (the repo is never touched).
import importlib.util,os,shutil,sys,tempfile,io,contextlib
ROOT='/home/claude/lsqvce'
spec=importlib.util.spec_from_file_location('audit_tree',os.path.join(ROOT,'dev/docs_src/audit_tree.py'))
A=importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
def run(check):
    A.FINDINGS.clear(); buf=io.StringIO()
    with contextlib.redirect_stdout(buf): check()
    return list(A.FINDINGS), buf.getvalue()
def case(name,check,mutate):
    d=tempfile.mkdtemp(prefix='swe_tests_inj_')
    shutil.copytree(ROOT,d,dirs_exist_ok=True,ignore=shutil.ignore_patterns('.git','__pycache__'))
    cwd=os.getcwd(); os.chdir(d)
    try:
        h,_=run(check); mutate(d); f,out=run(check)
    finally:
        os.chdir(cwd); shutil.rmtree(d,ignore_errors=True)
    print('%-62s healthy:%d after-fault:%d %s'%(name,len(h),len(f),'FIRES' if f else 'SILENT (blind)'))
def app(p,t):
    open(p,'a',encoding='utf-8').write(t)
T='dev/tests/test_signalCore.m'
case('10 duplicate local function in a TEST file (load error)',A.check_duplicate_functions,
     lambda d: app(os.path.join(d,T),'\nfunction makeCfg()\nend\nfunction makeCfg()\nend\n'))
case('10 duplicate local function in a PACKAGE file (control)',A.check_duplicate_functions,
     lambda d: app(os.path.join(d,'+lsqvce/+util/getfielddef.m'),'\nfunction zz()\nend\nfunction zz()\nend\n'))
case('7  unbalanced end in a TEST file',A.check_balance,
     lambda d: app(os.path.join(d,T),'\nfunction zzUnbalanced(tc)\n if true\nend\n'))
case('2  removed name used in TEST code',A.check_removed_names,
     lambda d: app(os.path.join(d,T),'\nfunction zz(tc)\n  x = coalesceOpt(1);\nend\n'))
case('3  nc as leading arg to adjustBlock in TEST code',A.check_signatures,
     lambda d: app(os.path.join(d,T),'\nfunction zz(tc)\n  adjustBlock(6, 1:3);\nend\n'))
case('15 wrong H1 line in a TEST file',A.check_h1_names,
     lambda d: open(os.path.join(d,T),'w').write(open(os.path.join(d,T)).read().replace('%TEST_SIGNALCORE','%TEST_WRONGNAME',1)))
case('15 wrong H1 line in a PACKAGE file (control)',A.check_h1_names,
     lambda d: open(os.path.join(d,'+lsqvce/+util/getfielddef.m'),'w').write(open(os.path.join(d,'+lsqvce/+util/getfielddef.m')).read().replace('%GETFIELDDEF','%WRONGNAME',1)))
case('9  spec field list (reads +lsqvce/adjustBlock.m, deleted)',A.check_spec_fields,lambda d: None)
case('14 id cited in a help block, raised only by a TEST file',A.check_error_ids,
     lambda d: (app(os.path.join(d,'+lsqvce/+util/getfielddef.m'),'\n%   lsqvce:zz:ghost  never raised by the package\n'),
                app(os.path.join(d,T),"\nfunction zz(tc)\n error('lsqvce:zz:ghost','x');\nend\n")))
case('16 tol: comment saying nothing (any text after tol:)',A.check_test_categories,
     lambda d: app(os.path.join(d,T),"\nfunction testZz(tc)\n% categories: precision\n    % tol:\n    tc.verifyEqual(1,1,'AbsTol',0.5);\nend\n"))
case('16 unit credited by a COMMENT mention only',A.check_test_categories,
     lambda d: None)
case('17 counts claim vs tree after adding a commented-out (%{ %}) test',A.check_readme_sync,
     lambda d: app(os.path.join(d,T),"\n%{\nfunction testHidden(tc)\nend\n%}\n"))
case('1  a test file under dev/tests/accum (stale path)',A.check_test_registration,lambda d: None)
