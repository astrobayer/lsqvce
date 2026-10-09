# Extension of audit check 27 to GROUP FIELDS: is each Options group field read anywhere outside its
# declaration (the defaults struct literal / help) -- engineArgs, fit, or the engine?
import re, glob, os
os.chdir('/home/claude/lsqvce')
src = open('+lsqvce/Options.m').read()
code = '\n'.join(l.split('%')[0] for l in src.split('\n'))
groups = {'Noise': ['Grid','Whitener','MeanType','Smooth','FreezeRedundancy','SkipThreshold'],
          'RobustOpts': ['Tuning','MaxIter','Tol'],
          'Diagnostics': ['Alpha','Bonferroni','Beta','RMin','RWeak','Test','RobustCheck','Store','NWorst'],
          'Parallel': ['Mode','DistributeDays','Threshold'], 'GaussNewton': ['TolF','Lambda0']}
pkg = ''
for f in glob.glob('+lsqvce/**/*.m', recursive=True):
    if f.endswith('+lsqvce/Options.m'): continue
    pkg += '\n'.join(l.split('%')[0] for l in open(f).read().split('\n'))
for g, fs in groups.items():
    for fld in fs:
        # a READ is <var>.<fld> (N.Grid, D.Alpha, obj.RobustOpts.Tuning, opts.Noise.Whitener)
        rd_opt = len(re.findall(r'\b\w+\.%s\b' % fld, code))
        rd_pkg = len(re.findall(r'(opts?|obj)\.%s\.%s\b' % (g, fld), pkg))
        print('%-12s %-17s reads in Options.m code: %d   opts.%s.%s reads in package: %d %s' %
              (g, fld, rd_opt, g, fld, rd_pkg, '<-- NO CONSUMER' if rd_opt + rd_pkg == 0 else ''))
