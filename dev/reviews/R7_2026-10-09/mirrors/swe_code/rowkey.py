# Mirror of lsqvce.util.rowWeightKey (asText: sprintf('%g', double(v))) vs makeLoader>dayKey ('%d')
def g(v): return '%g' % float(v)
def rowWeightKey(sys, day, ot, arc): return '%s|%s|%s|%s' % (g(sys), g(day), ot, g(arc))
for name, days in [('MJD 59396..59425', range(59396, 59426)),
                   ('datenum 739000..739029', range(739000, 739030)),
                   ('YYYYMMDD Sep 2026', range(20260901, 20260931)),
                   ('GPS-week*10+dow 2380x', range(23800, 23830)),
                   ('day-of-campaign 1e6+k', range(1000001, 1000031))]:
    keys = [rowWeightKey(1, d, 'kbr', 1) for d in days]
    dk = ['%d|kbr' % d for d in days]
    print('%-26s days=%d rowWeightKeys unique=%d  dayKey unique=%d  e.g. %s' %
          (name, len(keys), len(set(keys)), len(set(dk)), keys[0]))
print(rowWeightKey(1234567, 1, 'kbr', 1), rowWeightKey(1234568, 1, 'kbr', 1))
