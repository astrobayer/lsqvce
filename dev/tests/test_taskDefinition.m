function tests = test_taskDefinition
%TEST_TASKDEFINITION  Tests for the v7 declaration: TaskDefinition and its companions.
%
%     lsqvce.TaskDefinition    flat add-calls -> the engines' declaration struct
%     lsqvce.SystemDefinition  locals, streams (Build), priors of one system (Prepare)
%     lsqvce.Noise             one constructor per temporal law
%     lsqvce.DayContext        what Prepare and Build know about their day
%     lsqvce.DayData           the arcs a Build returns
%
%   Claims under test: (1) a task compiles to exactly the struct a user
%   would write by hand, and its fit equals that struct's fit BITWISE --
%   one solve and a noise run (the same arithmetic in the same order, so no
%   tolerance); (2) local partials by name (Local=) and as a matrix (B=)
%   give the bitwise same fit; names at once equal names one at a time;
%   (3) every refusal fires its documented identifier, at the call that
%   makes the mistake or while the day is read; (4) the context announces
%   the day, the schedule, the columns and the tables as the engines' layout
%   has them, including a data gap and a validity-limited local; (5) the
%   generated wrappers cost little next to the fit (a ratio budget).
%
%   Run: runtests('dev/tests/test_taskDefinition') with the toolbox root on
%   the path, or through ci_gate.
%
%   Created 2026-11-02 (CET, UTC+01:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- level v7, WP-160.
    tests = functiontests(localfunctions);
end

%==========================================================================
function testCompiledDeclarationEqualsTheHandWrittenStruct(tc)
% categories: contract, reference
    % The task of the fixture compiles to the struct a v6 user writes by hand,
    % field by field; only the two handles differ (generated wrappers).
    F = fixture();
    tc.verifyInstanceOf(F.task, ?lsqvce.TaskDefinition);
    tc.verifyInstanceOf(F.sys, ?lsqvce.SystemDefinition);
    tc.verifyInstanceOf(F.noiseA, ?lsqvce.Noise);
    dT = F.task.declaration();
    dH = F.data;
    tc.verifyEqual(dT.commonParams, dH.commonParams, 'the common template table');
    tc.verifyEqual(dT.days, dH.days, 'the calendar');
    tc.verifyEqual(dT.paramColumns, dH.paramColumns, 'validity and name columns');
    tc.verifyEqual(dT.priors, dH.priors, 'the common priors (On -> group, L)');
    sT = rmfield(dT.systems, {'buildDay', 'prepareDay'});
    sH = rmfield(dH.systems, {'buildDay', 'prepareDay'});
    tc.verifyEqual(sT, sH, 'the system: covariance, localParams, arcs, priors, name');
    tc.verifyEqual(fieldnames(dT.systems.covariance), {'a'; 'p'}, 'streams in the order added');
    tc.verifyClass(dT.systems.buildDay, 'function_handle');
    tc.verifyClass(dT.systems.prepareDay, 'function_handle');
end

function testCompiledTaskFitsBitwiseLikeTheStruct(tc)
% categories: reference, precision
    % THE equivalence claim of the level: the same declaration through the
    % generated wrappers gives the same numbers BITWISE -- one solve and a
    % noise run (levels estimated, three passes).
    % tol: exact (isequal) -- the same arithmetic in the same order
    F = fixture();
    tc.verifyInstanceOf(F.sys, ?lsqvce.SystemDefinition);
    tc.verifyInstanceOf(F.noiseA, ?lsqvce.Noise);
    dT = F.task.declaration();
    mcC = ?lsqvce.DayContext;
    mcD = ?lsqvce.DayData;
    tc.verifySubstring(func2str(dT.systems.buildDay), [mcC.Name '.runBuild'], ...
        'the streams compile to the DayContext wrapper');
    tc.verifyEqual(mcD.Name, 'lsqvce.DayData');
    o1 = lsqvce.Options(Estimate="none", Verbose=false);
    R1 = lsqvce.fit(F.task, o1);
    H1 = lsqvce.engine.fitDeclaration(F.data, o1);
    tc.verifyTrue(isequal(R1.x, H1.x), 'one solve: x bitwise');
    tc.verifyTrue(isequal(R1.sigma0, H1.sigma0), 'one solve: sigma0 bitwise');
    tc.verifyTrue(isequal(R1.sol.Cxx, H1.sol.Cxx), 'one solve: Cxx bitwise');
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');      % MaxPasses 3 on purpose
    o2 = lsqvce.Options(Estimate="noise", MaxPasses=3, Seed=7, Verbose=false);
    R2 = lsqvce.fit(F.task, o2);
    H2 = lsqvce.engine.fitDeclaration(F.data, o2);
    tc.verifyTrue(isequal(R2.x, H2.x), 'noise run: x bitwise');
    tc.verifyTrue(isequal(R2.sigma0, H2.sigma0), 'noise run: sigma0 bitwise');
    tc.verifyTrue(isequal(R2.sol.Cxx, H2.sol.Cxx), 'noise run: Cxx bitwise');
    tc.verifyEqual(R2.stopReason, H2.stopReason);
    % and the fit is a fit: the planted commons within 4 formal sigma
    se = sqrt(diag(R1.sol.Cxx));
    % tol: 4 formal sigma -- 600 rows of white noise, a statistical band
    tc.verifyLessThan(max(abs(R1.x(1:2) - [0.3; -0.8]) ./ se(1:2)), 4);
end

function testLocalByNameEqualsLocalAsMatrix(tc)
% categories: reference, precision
    % D4: the fixture with Local={names, P} and its twin with B= built from
    % ctx.cols give the bitwise same fit; and one Local row of several
    % names equals one row per name (the same B, the same fit).
    % tol: exact (isequal) -- the same B, the same arithmetic
    byName = fixture("Local");
    byB = fixture("B");
    mc = ?lsqvce.DayData;
    tc.verifyEqual(mc.Name, 'lsqvce.DayData');
    mcC = ?lsqvce.DayContext;
    tc.verifyEqual(mcC.Name, 'lsqvce.DayContext');
    tc.verifyInstanceOf(byB.task, ?lsqvce.TaskDefinition);
    o = lsqvce.Options(Estimate="none", Verbose=false);
    R1 = lsqvce.fit(byName.task, o);
    R2 = lsqvce.fit(byB.task, o);
    tc.verifyTrue(isequal(R1.x, R2.x) && isequal(R1.sol.Cxx, R2.sol.Cxx), ...
        'Local= and B= give the same solution bitwise');
    T1 = pairTask("together");
    T2 = pairTask("apart");
    tc.verifyTrue(isequal(lsqvce.fit(T1, o).x, lsqvce.fit(T2, o).x), ...
        'names in one Local row = one row per name');
end

function testNamesAtOnceEqualNamesOneAtATime(tc)
% categories: vectorisation
    % BITWISE: addCommon / addLocal of a name list compile to the same
    % tables as one call per name.
    t1 = pairTask("together");
    t2 = lsqvce.TaskDefinition(Days=1:2);
    t2.addCommon("u");  t2.addCommon("v");
    s2 = t2.addSystem("S", Arcs=1);
    s2.addLocal("i", Scope="arc");  s2.addLocal("j", Scope="arc");
    s2.addStream("y", Build=@(d, ctx) pairBuild(d, ctx, "apart"), Noise=lsqvce.Noise.white(1));
    tc.verifyInstanceOf(s2, ?lsqvce.SystemDefinition);
    tc.verifyInstanceOf(t2, ?lsqvce.TaskDefinition);
    d1 = t1.declaration();  d2 = t2.declaration();
    tc.verifyEqual(d1.commonParams, d2.commonParams);
    tc.verifyEqual(d1.systems.localParams, d2.systems.localParams);
    tc.verifyEqual(d1.systems.covariance, d2.systems.covariance);
end

function testInMemoryTasksFitLikeTheirStructs(tc)
% categories: contract, reference
    % D6: fromMatrices / fromBatches compile to the matrices / batched kinds;
    % a prior goes by L. The fits are bitwise those of the structs.
    rs = RandStream('mt19937ar', 'Seed', 3);
    A = randn(rs, 50, 3);  l = A * [1; 2; 3] + 0.1 * randn(rs, 50, 1);  p = 1 + rand(rs, 50, 1);
    t = lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p);
    t.addPrior("tie", L=[1 -1 0], Sigma=0.5);
    tc.verifyInstanceOf(t, ?lsqvce.TaskDefinition);
    d = struct('A', A, 'l', l, 'p', p, 'priors', struct('group', [], 'L', [1 -1 0], 'x0', 0, ...
               'sigma', 0.5, 'from', 1, 'name', 'tie', 'estimate', "none"));
    tc.verifyEqual(t.declaration(), d);
    o = lsqvce.Options(Verbose=false);
    % tol: exact (isequal) -- the same struct
    tc.verifyTrue(isequal(lsqvce.fit(t, o).x, lsqvce.engine.fitDeclaration(d, o).x));
    b = lsqvce.TaskDefinition.fromBatches(@(k) deal(A(k:5:end, :), l(k:5:end), p(k:5:end)), ...
                                          NumBatches=5, NumUnknowns=3);
    tc.verifyEqual(b.declaration(), struct('batchFcn', b.declaration().batchFcn, 'numBatches', 5, ...
                                           'numUnknowns', 3));
    % tol: 1e-12 relative -- batched against one block: cond(A) ~ 2, rounding of the fold order
    tc.verifyEqual(lsqvce.fit(b, o).x, lsqvce.engine.fitDeclaration(struct('A', A, 'l', l, 'p', p), o).x, 'RelTol', 1e-12);
    g = lsqvce.TaskDefinition.fromBatches(@(k, x) deal(A(k:5:end, :), l(k:5:end) - A(k:5:end, :) * x, ...
                                                       p(k:5:end)), NumBatches=5, X0=zeros(3, 1));
    tc.verifyEqual(g.declaration().x0, zeros(3, 1), 'X0 makes the nonlinear form');
    tc.verifyError(@() t.addCommon("c"), 'lsqvce:TaskDefinition:wrongKind');
    tc.verifyError(@() t.addSystem("S"), 'lsqvce:TaskDefinition:wrongKind');
    tc.verifyError(@() t.addPrior("p", On="^c", Sigma=1), 'lsqvce:TaskDefinition:wrongKind');
    tc.verifyError(@() t.layout(), 'lsqvce:TaskDefinition:wrongKind');
    tc.verifyError(@() lsqvce.TaskDefinition.fromBatches(@(k) k), 'lsqvce:TaskDefinition:missingArgument');
    shown = evalc('disp(t)');
    echoWarnings(shown);   % WP-145: a capture passes its warnings on
    tc.verifySubstring(shown, 'matrices');
end

function testDeclarationRefusesWhatItCannotCompile(tc)
% categories: contract
    % Every declaration-time identifier fires at the call that makes the
    % mistake (or at declaration(), for what only the whole task can know).
    t = lsqvce.TaskDefinition(Days=1:2);
    t.addCommon(["a" "b"]);
    tc.verifyError(@() t.addCommon("a"), 'lsqvce:TaskDefinition:duplicateName');
    tc.verifyError(@() t.addCommon(["c" "c"]), 'lsqvce:TaskDefinition:duplicateName');
    tc.verifyError(@() t.addCommon("c", Scope="arc"), 'MATLAB:validators:mustBeMember');
    tc.verifyError(@() t.addCommon("c", Valid=[3 2]), 'lsqvce:TaskDefinition:badValid');
    tc.verifyError(@() t.addCommon("c", Valid=[1 2 3]), 'lsqvce:TaskDefinition:badValid');
    tc.verifyError(@() t.addCommon(["c" "d"], Info=table(1)), 'lsqvce:TaskDefinition:infoRows');
    tc.verifyError(@() t.addSystem("S", Prepare=42), 'lsqvce:TaskDefinition:badHandle');
    tc.verifyError(@() t.addSystem("S", Prepare=@(d) d), 'lsqvce:TaskDefinition:badHandle');
    s = t.addSystem("S");
    tc.verifyInstanceOf(s, ?lsqvce.SystemDefinition);
    tc.verifyError(@() t.addSystem("S"), 'lsqvce:TaskDefinition:duplicateName');
    s.addLocal("i", Scope="day");
    tc.verifyError(@() s.addLocal("i"), 'lsqvce:TaskDefinition:duplicateName');
    tc.verifyError(@() s.addStream("bad name", Build=@(d, c) c, Noise=lsqvce.Noise.white(1)), ...
                   'lsqvce:TaskDefinition:badName');
    tc.verifyError(@() s.addStream("y", Noise=lsqvce.Noise.white(1)), 'lsqvce:TaskDefinition:missingArgument');
    tc.verifyError(@() s.addStream("y", Build=@(d, c) c), 'lsqvce:TaskDefinition:missingArgument');
    tc.verifyError(@() s.addStream("y", Build=@(d) d, Noise=lsqvce.Noise.white(1)), ...
                   'lsqvce:TaskDefinition:badHandle');
    tc.verifyError(@() s.addStream("y", Build=@(d, c) c, Noise=lsqvce.Noise.white(1), Estimate="all"), ...
                   'MATLAB:validators:mustBeMember');
    tc.verifyError(@() t.declaration(), 'lsqvce:TaskDefinition:incomplete');    % a system without a stream
    s.addStream("y", Build=@(d, c) c, Noise=lsqvce.Noise.white(1));
    tc.verifyError(@() s.addStream("y", Build=@(d, c) c, Noise=lsqvce.Noise.white(1)), ...
                   'lsqvce:TaskDefinition:duplicateName');
    tc.verifyError(@() t.addPrior("p", Sigma=1), 'lsqvce:TaskDefinition:priorTarget');
    tc.verifyError(@() t.addPrior("p", On="^a", L=[1 0], Sigma=1), 'lsqvce:TaskDefinition:priorTarget');
    tc.verifyError(@() t.addPrior("p", On="^a"), 'lsqvce:TaskDefinition:missingArgument');
    t.addPrior("p", On="^a", Sigma=1);
    tc.verifyError(@() t.addPrior("p", On="^b", Sigma=1), 'lsqvce:TaskDefinition:duplicateName');
    s.addLocal("a", Scope="day");                                               % named like a common
    tc.verifyError(@() t.declaration(), 'lsqvce:TaskDefinition:duplicateName');
    % arcs on some streams only; no calendar; nothing declared
    u = lsqvce.TaskDefinition(Days=1);
    u.addCommon("c");
    w = u.addSystem("W");
    w.addStream("y", Build=@(d, c) c, Noise=lsqvce.Noise.white(1), Arcs=1);
    w.addStream("z", Build=@(d, c) c, Noise=lsqvce.Noise.white(1));
    tc.verifyError(@() u.declaration(), 'lsqvce:TaskDefinition:badArcs');
    v = lsqvce.TaskDefinition();
    v.addCommon("c");
    v.addSystem("V").addStream("y", Build=@(d, c) c, Noise=lsqvce.Noise.white(1));
    tc.verifyError(@() v.declaration(), 'lsqvce:TaskDefinition:incomplete');   % no calendar
    tc.verifyError(@() lsqvce.TaskDefinition(Days=1).declaration(), 'lsqvce:TaskDefinition:incomplete');
    tc.verifyError(@() lsqvce.fit(42), 'MATLAB:validators:mustBeA');
    % WP-162: the struct left the user tier -- refused by name, with the way out
    d = struct('A', eye(2), 'l', [1; 2]);
    tc.verifyError(@() lsqvce.fit(d), 'lsqvce:fit:structDeclaration');
    tc.verifyError(@() lsqvce.layout(d), 'lsqvce:layout:structDeclaration');
end

function testDayContextRefusesWhatItCannotPlace(tc)
% categories: contract, robustness
    % Every read-time identifier fires while the day is read, unchanged
    % through fit (an lsqvce:* error passes the engines' wrapper as it is).
    mc = ?lsqvce.DayContext;
    tc.verifyEqual(mc.Name, 'lsqvce.DayContext');
    mcD = ?lsqvce.DayData;
    tc.verifyEqual(mcD.Name, 'lsqvce.DayData');
    o = lsqvce.Options(Verbose=false);
    cases = {"unknownName", 'lsqvce:DayContext:unknownName'; ...
             "unknownCols", 'lsqvce:DayContext:unknownName'; ...
             "arcRequired", 'lsqvce:DayContext:arcRequired'; ...
             "arcNotScheduled", 'lsqvce:DayContext:arcRequired'; ...
             "mixed", 'lsqvce:DayContext:mixedTables'; ...
             "badReturn", 'lsqvce:DayContext:badReturn'; ...
             "localAndB", 'lsqvce:DayContext:localAndB'; ...
             "localOneColumn", 'lsqvce:DayContext:badLocal'; ...
             "localWrongSize", 'lsqvce:DayContext:badLocal'; ...
             "localTwice", 'lsqvce:DayContext:badLocal'; ...
             "localNotText", 'lsqvce:DayContext:badLocal'; ...
             "arcUnscheduled", 'lsqvce:DayContext:badArc'; ...
             "arcTwice", 'lsqvce:DayContext:badArc'};
    for k = 1:size(cases, 1)
        tc.verifyError(@() lsqvce.fit(tinyTask(cases{k, 1}, true), o), cases{k, 2}, ...
                       sprintf('case %s', cases{k, 1}));
    end
    tc.verifyError(@() lsqvce.fit(tinyTask("arcWithoutSchedule", false), o), 'lsqvce:DayContext:badArc');
    tc.verifyError(@() lsqvce.fit(tinyTask("dayNotOnDay", true), o), 'lsqvce:DayContext:unknownName');
    R = lsqvce.fit(tinyTask("ok", true), o);
    tc.verifyTrue(all(isfinite(R.x)), 'the good twin of every case fits');
end

function testContextAnnouncesTheDay(tc)
% categories: contract, robustness
    % The context of each call: the day, the names, Settings without the
    % reserved field, Prepared, the counts, the stream's schedule (a stream
    % resolver that differs from the system's, and a gap), the columns as
    % the engines' layout has them, the tables with OnDay (validity).
    F = fixture("Local", true);
    tc.verifyInstanceOf(F.task, ?lsqvce.TaskDefinition);
    tc.verifyInstanceOf(F.sys, ?lsqvce.SystemDefinition);
    R = lsqvce.fit(F.task, lsqvce.Options(Verbose=false));
    tc.verifyTrue(all(isfinite(R.x)), 'the fit runs with a gap in stream q');
    seen = F.recorder;
    P = seen('prepare|2');
    tc.verifyEqual(P.Stream, "", 'Prepare has no stream');
    tc.verifyEmpty(P.Prepared, 'Prepare sees no Prepared');
    tc.verifyEqual(P.Settings.note, 'mine', 'Prepare sees the Settings');
    tc.verifyEqual(P.ArcIds, [1 2 3], 'Prepare sees the system''s arcs of the day (a, p: 1 2; q: 1 3)');
    B = seen('a|2');
    tc.verifyEqual(B.Day, 2);
    tc.verifyEqual([B.System, B.Stream], ["S", "a"]);
    tc.verifyFalse(isfield(B.Settings, 'prepared'), 'Settings without the reserved field');
    tc.verifyEqual(B.Prepared.day, 2, 'Prepared is this day''s');
    lay = F.task.layout();
    tc.verifyEqual(B.NumCommon, lay.nc);
    tc.verifyEqual(B.NumLocal, lay.systems(1).nl);
    tc.verifyEqual(B.ArcIds, [1 2], 'the system schedule');
    rowIc = find([lay.systems(1).local.name] == "ic", 1);
    rowIc = lay.systems(1).local(rowIc).row;
    tc.verifyEqual(B.colsIc2, lay.systems(1).col(rowIc, 2, 2), 'ctx.cols = the layout, arc 2');
    jGd = find(arrayfun(@(c) c.name == "gd" && isequal(c.day, 2), lay.common), 1);
    rowGd = lay.common(jGd).row;
    tc.verifyEqual(B.colsGd, lay.colCommon(rowGd, 2), 'a day common moves with the day');
    tc.verifyEqual(B.localOn, struct('ic', true, 'bias', true, 'cal', true, 'man', true), ...
                   'man (Valid [2.25 2.75]) is on day 2');
    B1 = seen('a|1');
    tc.verifyFalse(B1.localOn.man, 'and not on day 1');
    Q = seen('q|2');
    tc.verifyEqual(Q.ArcIds, [1 3], 'the stream''s own resolver');
    tc.verifyEqual(lsqvce.layout(F.task).nc, lay.nc, 'lsqvce.layout takes a task');
    out = evalc('disp(F.task)');
    echoWarnings(out);     % WP-145: a capture passes its warnings on
    tc.verifySubstring(out, 'system "S"');
    tc.verifySubstring(out, 'prior');
end

function testToTableOfATaskIsOneTable(tc)
% categories: contract
    % WP-161: R.toTable(task) is ONE table, a row per instantiated parameter
    % (commons in A's order, then each system's locals in B's order), with
    % names, scopes, day / arc of the instance and the estimate with its
    % standard deviation -- the same numbers as R.x and the cell tables.
    F = fixture();
    tc.verifyInstanceOf(F.task, ?lsqvce.TaskDefinition);
    R = lsqvce.fit(F.task, lsqvce.Options(Verbose=false));
    T = R.toTable(F.task);
    % WP-162: the struct's write-back is the engine tier's (hidden declarationTables)
    tc.verifyError(@() R.toTable(F.data), 'lsqvce:Results:structDeclaration');
    lay = F.task.layout();
    tc.verifyClass(T, 'table');
    tc.verifyEqual(T.Properties.VariableNames, {'System', 'Name', 'Scope', 'Day', 'Arc', 'Estimate', 'Std'});
    tc.verifyEqual(height(T), lay.nc + lay.systems(1).nl, 'one row per instantiated parameter');
    tc.verifyEqual(T.Estimate(1:lay.nc), R.x, 'the commons are R.x, in order');
    tc.verifyEqual(T.Std(1:lay.nc), sqrt(diag(R.sol.Cxx)), 'and their standard deviations');
    tc.verifyEqual(T.Name(1:lay.nc).', [lay.common.name], 'named as the layout names them');
    tc.verifyEqual(T.System(lay.nc + 1), "S");
    gd = T(T.Name == "gd", :);
    tc.verifyEqual(gd.Day.', 1:3, 'a day common once per day');
    man = T(T.Name == "man", :);
    tc.verifyEqual(man.Day, 2, 'the validity-limited local on its one day');
    ic = T(T.Name == "ic", :);
    tc.verifyEqual(height(ic), 6, 'an arc local per (day, arc)');
    tc.verifyTrue(all(isnan(T.Arc(T.Scope ~= "arc"))), 'Arc is NaN off the arc scope');
end

function testTaskOverheadIsSmall(tc)
% categories: speed
    % The generated wrappers (a DayContext and a DayData per build, columns
    % by name) against the hand-written buildDay of the same day, one solve.
    F = fixture();
    tc.verifyInstanceOf(F.task, ?lsqvce.TaskDefinition);
    o = lsqvce.Options(Estimate="none", Verbose=false);
    rec = timePair(@() lsqvce.fit(F.task, o), @() lsqvce.engine.fitDeclaration(F.data, o), 'Repeats', 9, ...
                   'Label', 'task vs struct, one solve');
    tc.log(1, sprintf('WP-160 overhead: ratio %.3f, band %.3f-%.3f, tA %.4f s, tB %.4f s', ...
                      rec.ratio, rec.band, rec.tA, rec.tB));
    % tol: 2.0 -- measured 1.43-1.63 over four runs on R2026a (PCWIN64, WP-160):
    % ~1.7 ms per build of this 600-row fixture (ctx.Local's table 0.6 ms,
    % addArc and cols the rest); a real Build (files, orbit integration) costs
    % orders of magnitude more, so the ratio is the worst case
    tc.verifyLessThan(rec.ratio, 2.0, sprintf('ratio %.2f (band %.2f-%.2f)', rec.ratio, rec.band));
end

function testNoiseIsTheEnginesSpecification(tc)
% categories: contract, reference, robustness
    % Each law equals the engines' covariance specification field by field;
    % the validators refuse what is not a law; a per-row sigma is a fixed
    % structure, so the stream's Estimate="law" is refused by the engines.
    rho = @(dt, ~) exp(-abs(dt) / 30);
    cf = @(t, ~) eye(numel(t));
    tc.verifyEqual(lsqvce.Noise.white(2).spec(), struct('type', 'diag', 'sigma', 2));
    tc.verifyEqual(lsqvce.Noise.arma([1 -0.5], 1.2, MA=[1 0.3]).spec(), ...
                   struct('type', 'arma', 'ar', [1 -0.5], 'ma', [1 0.3], 'sigma', 1.2));
    tc.verifyEqual(lsqvce.Noise.psd([0 1].', [2 1].', MaxDenseArc=500).spec(), ...
                   struct('type', 'psd', 'f', [0; 1], 'S', [2; 1], 'maxDenseArc', 500));
    tc.verifyEqual(lsqvce.Noise.acf(rho, 3).spec(), struct('type', 'acf', 'acf', rho, 'sigma', 3));
    tc.verifyEqual(lsqvce.Noise.var(0.5 * eye(2), eye(2)).spec(), ...
                   struct('type', 'var', 'ar', 0.5 * eye(2), 'sigma', eye(2)));
    tc.verifyEqual(lsqvce.Noise.custom(cf).spec(), struct('type', 'custom', 'fcn', cf));
    tc.verifyEqual(lsqvce.Noise.white(2).Law, "white");
    tc.verifyError(@() lsqvce.Noise.white(-1), 'MATLAB:validators:mustBePositive');
    tc.verifyError(@() lsqvce.Noise.white(Inf), 'MATLAB:validators:mustBeFinite');
    tc.verifyError(@() lsqvce.Noise.acf(3, 1), 'MATLAB:validation:UnableToConvert');
    tc.verifyError(@() lsqvce.Noise.white(1, MaxDenseArc=1.5), 'MATLAB:validators:mustBeInteger');
    % a per-row sigma vector: a level, no shape
    t = lsqvce.TaskDefinition(Days=1);
    t.addCommon("c");
    t.addSystem("S").addStream("y", Build=@(d, ctx) pairBuild(d, ctx, "plain"), ...
        Noise=lsqvce.Noise.white(ones(20, 1)), Estimate="law");
    tc.verifyError(@() lsqvce.fit(t, Estimate="noise", Verbose=false), 'lsqvce:fit:lawOnFixedStructure');
end

%==========================================================================
% fixtures
%==========================================================================
function F = fixture(localMode, record)
%FIXTURE  One system "S" over days 1:3: commons g1, g2 ('common') and gd
%   ('day'); locals ic ('arc'), bias ('day'), cal ('system'), man ('day',
%   Valid [2.25 2.75] -> day 2 only); a Prepare (makeDay); stream a (scalar,
%   white 0.5, one excluded pair of rows on day 1) and stream p (3 axes, a
%   per-axis Sigma); priors On "^g" and L on the day commons (common), On
%   "^cal" (local). F.task is the TaskDefinition, F.data the struct a v6 user
%   writes by hand for the same problem. LOCALMODE "Local" (default) gives
%   the partials by name, "B" as a matrix from ctx.cols. RECORD adds a
%   recorder (F.recorder) and a third stream q with its own resolver.
    if nargin < 1, localMode = "Local"; end
    if nargin < 2, record = false; end
    rec = containers.Map('KeyType', 'char', 'ValueType', 'any');
    task = lsqvce.TaskDefinition(Days=1:3, Name="fixture");
    task.addCommon(["g1" "g2"]);
    task.addCommon("gd", Scope="day");
    task.addPrior("loose", On="^g", Sigma=10);
    task.addPrior("dayTie", L=[zeros(3, 2), eye(3)], Value=0.2, Sigma=1);
    settings = struct();
    if record, settings.note = 'mine'; end
    sys = task.addSystem("S", Prepare=@(day, ctx) recordPrepare(day, ctx, rec), Settings=settings, Arcs=2);
    sys.addLocal("ic", Scope="arc");
    sys.addLocal("bias", Scope="day");
    sys.addLocal("cal", Scope="system");
    sys.addLocal("man", Scope="day", Valid=[2.25 2.75]);
    sys.addPrior("cal", On="^cal", Value=1, Sigma=0.1);
    noiseA = lsqvce.Noise.white(0.5);
    sys.addStream("a", Build=@(day, ctx) taskBuild(day, ctx, 'a', localMode, rec), Noise=noiseA, ...
                  Estimate="scale", ArcFactors=false);
    sys.addStream("p", Build=@(day, ctx) taskBuild(day, ctx, 'p', localMode, rec), ...
                  Noise=lsqvce.Noise.white(1), Estimate="scale", Axes=3);
    if record
        sys.addStream("q", Build=@(day, ctx) qBuild(day, ctx, rec), Noise=lsqvce.Noise.white(0.5), ...
                      Arcs=@(day, ~) qSchedule(day));
    end
    F = struct('task', task, 'sys', sys, 'noiseA', noiseA, 'data', handStruct(), 'recorder', rec);
end

function data = handStruct()
%HANDSTRUCT  The v6 struct of the fixture, written as a user writes it.
    data = struct();
    data.commonParams = {'common', 'g1', [], []; 'common', 'g2', [], []; 'day', 'gd', [], []};
    data.days = 1:3;
    s = struct();
    s.buildDay = @(day, obsType, common, local, procOpt, info) handBuild(obsType, procOpt, info);
    s.covariance = struct('a', struct('type', 'diag', 'sigma', 0.5, 'estimate', "scale", 'arcFactors', false), ...
                          'p', struct('type', 'diag', 'sigma', 1, 'estimate', "scale"));
    s.localParams = {'arc', 'ic', [], []; 'day', 'bias', [], []; 'system', 'cal', [], []; ...
                     'day', 'man', 2.25, 2.75};
    s.arcs = 2;
    s.days = [];
    s.id = [];
    s.name = 'S';
    s.processingOptions = struct();
    s.priors = struct('group', '^cal', 'L', [], 'x0', 1, 'sigma', 0.1, 'from', 1, 'name', 'cal', ...
                      'estimate', "none");
    s.prepareDay = @(day, common, local, procOpt, info) makeDay(day);
    data.systems = s;
    data.paramColumns = struct('name', 2, 'start', 3, 'end', 4);
    data.priors = struct('group', {'^g', []}, 'L', {[], [zeros(3, 2), eye(3)]}, 'x0', {0, 0.2}, ...
                         'sigma', {10, 1}, 'from', {1, 1}, 'name', {'loose', 'dayTie'}, ...
                         'estimate', {"none", "none"});
end

function P = makeDay(day)
%MAKEDAY  The content of one day, independent of any layout: per stream and
%   arc the common partials (g1 g2 | gd), the local partials, o - c, epochs,
%   per-axis sigmas (p) and excluded rows (a, day 1).
    rs = RandStream('mt19937ar', 'Seed', 100 + day);
    gd = [0.1 0.2 0.3];
    P = struct('day', day);
    for w = 'ap'
        if w == 'a', m = 40; else, m = 60; end
        for a = 1:2
            q = struct();
            q.Ac = randn(rs, m, 2);
            q.Ad = randn(rs, m, 1);
            q.ic = 1 + (0:m-1).' / m;
            q.bias = ones(m, 1);
            q.cal = 0.5 * randn(rs, m, 1);
            q.man = double(day == 2) * cos((1:m).' / 7);
            q.t = (day - 1) * 86400 + (a - 1) * 43200 + 30 * (0:m-1).';
            if w == 'p'
                q.sig = 0.5 + rand(rs, m / 3, 3);
                e = reshape(q.sig.', [], 1) .* randn(rs, m, 1);
            else
                q.sig = [];
                e = 0.5 * randn(rs, m, 1);
            end
            q.l = q.Ac * [0.3; -0.8] + q.Ad * gd(day) + q.ic * (0.05 * a) + q.bias * (0.1 * day) ...
                  + q.cal * 1.05 + q.man * 0.7 + e;
            q.ex = false(m, 1);
            if w == 'a' && day == 1 && a == 1, q.ex(5:6) = true; end
            for f = {'Ac', 'Ad', 'ic', 'bias', 'cal', 'man', 'l'}
                q.(f{1})(q.ex, :) = 0;
            end
            P.(w)(a) = q;
        end
    end
end

function P = recordPrepare(day, ctx, rec)
%RECORDPREPARE  The fixture's Prepare: makeDay, and what the context said.
    rec(sprintf('prepare|%d', day)) = struct('Stream', ctx.Stream, 'Prepared', ctx.Prepared, ...
                                             'Settings', ctx.Settings, 'ArcIds', ctx.ArcIds);
    P = makeDay(day);
end

function out = taskBuild(day, ctx, w, localMode, rec)
%TASKBUILD  The fixture's Build in the v7 way: columns by name.
    P = ctx.Prepared;
    out = ctx.newDay();
    L = ctx.Local;
    onMan = L.OnDay(L.Name == "man");
    for a = 1:2
        q = P.(w)(a);
        m = numel(q.l);
        A = zeros(m, ctx.NumCommon);
        A(:, ctx.cols(["g1" "g2"])) = q.Ac;
        A(:, ctx.cols("gd")) = q.Ad;
        loc = {"ic", q.ic; "bias", q.bias; "cal", q.cal};
        if onMan, loc(end + 1, :) = {"man", q.man}; end %#ok<AGROW>
        ex = [];
        if w == 'a', ex = q.ex; end
        if localMode == "B"
            B = sparse(m, ctx.NumLocal);
            for r = 1:size(loc, 1)
                B(:, ctx.cols(loc{r, 1}, Arc=a)) = loc{r, 2};
            end
            out.addArc(A, q.l, q.t, B=B, Sigma=q.sig, Excluded=ex);
        else
            out.addArc(A, q.l, q.t, Local=loc, Sigma=q.sig, Excluded=ex);
        end
    end
    if rec.Count >= 0
        on = struct();
        for r = 1:height(L), on.(char(L.Name(r))) = L.OnDay(r); end
        rec(sprintf('%s|%d', w, day)) = struct('Day', ctx.Day, 'System', ctx.System, ...
            'Stream', ctx.Stream, 'Settings', ctx.Settings, 'Prepared', P, 'NumCommon', ctx.NumCommon, ...
            'NumLocal', ctx.NumLocal, 'ArcIds', ctx.ArcIds, 'colsIc2', ctx.cols("ic", Arc=2), ...
            'colsGd', ctx.cols("gd"), 'localOn', on);
    end
end

function D = handBuild(obsType, procOpt, info)
%HANDBUILD  The same day as a v6 user writes it: offsets from info.
    P = procOpt.prepared;
    w = obsType;
    [A, B, l, t, C, E] = deal(cell(1, 2));
    for a = 1:2
        q = P.(w)(a);
        m = numel(q.l);
        Aa = zeros(m, info.nc);
        Aa(:, 1:2) = q.Ac;
        Aa(:, info.commonDayOffset + find(info.commonDayRows == 3)) = q.Ad;
        Ba = sparse(m, info.nl);
        Ba(:, info.offsetOf(info.arcId(a)) + find(info.arcRows == 1)) = q.ic;
        Ba(:, info.dayOffset + find(info.dayRows == 2)) = q.bias;
        Ba(:, info.sysOffset + find(info.sysRows == 3)) = q.cal;
        if any(info.dayRows == 4)
            Ba(:, info.dayOffset + find(info.dayRows == 4)) = q.man;
        end
        A{a} = Aa;  B{a} = Ba;  l{a} = q.l;  t{a} = q.t;  C{a} = q.sig;  E{a} = q.ex;
    end
    D = struct('A', {A}, 'B', {B}, 'l', {l}, 't', {t});
    if w == 'p'
        D.Cep = C;
        D.nAxis = 3;
    else
        D.excluded = E;
    end
end

function ids = qSchedule(day)
%QSCHEDULE  Stream q's own schedule: arcs 1 and 3 on day 2, arc 1 otherwise.
    if day == 2, ids = [1 3]; else, ids = 1; end
end

function out = qBuild(day, ctx, rec)
%QBUILD  Stream q: on day 2 only arc 3 carries data (arc 1 is a gap).
    rec(sprintf('q|%d', day)) = struct('ArcIds', ctx.ArcIds);
    out = ctx.newDay();
    m = 30;
    rs = RandStream('mt19937ar', 'Seed', 500 + day);
    A = zeros(m, ctx.NumCommon);
    A(:, ctx.cols(["g1" "g2"])) = randn(rs, m, 2);
    ic = {"ic", 1 + (0:m-1).' / m};                    % the arc's own local is observed
    if day == 2
        out.addArc(A, A(:, 1:2) * [0.3; -0.8] + 0.5 * randn(rs, m, 1), 1e5 + 30 * (0:m-1).', ...
                   Local=ic, Arc=3);
    else
        out.addArc(A, A(:, 1:2) * [0.3; -0.8] + 0.5 * randn(rs, m, 1), 2e5 * day + 30 * (0:m-1).', ...
                   Local=ic);
    end
end

function task = pairTask(mode)
%PAIRTASK  Days 1:2, commons u v, arc locals i j, stream y through pairBuild(MODE).
    task = lsqvce.TaskDefinition(Days=1:2);
    task.addCommon(["u" "v"]);
    sys = task.addSystem("S", Arcs=1);
    sys.addLocal(["i" "j"], Scope="arc");
    sys.addStream("y", Build=@(d, ctx) pairBuild(d, ctx, mode), Noise=lsqvce.Noise.white(1));
end

function out = pairBuild(day, ctx, mode)
%PAIRBUILD  Two commons u, v and two arc locals i, j, given together or apart
%   ("together" / "apart"); "plain" is one arc of 20 rows without locals.
    rs = RandStream('mt19937ar', 'Seed', 900 + day);
    m = 20;
    out = ctx.newDay();
    A = randn(rs, m, ctx.NumCommon);
    l = A * (1:ctx.NumCommon).' + 0.1 * randn(rs, m, 1);
    t = 1e5 * day + 10 * (0:m-1).';
    switch mode
        case "plain"
            out.addArc(A, l, t);
        case "together"
            out.addArc(A, l, t, Local={["i" "j"], [ones(m, 1), (1:m).' / m]});
        otherwise
            out.addArc(A, l, t, Local={"i", ones(m, 1); "j", (1:m).' / m});
    end
end

function task = tinyTask(mode, scheduled)
%TINYTASK  One day, common c, locals ic ('arc', with a schedule) and b
%   ('day'), one stream whose Build makes the mistake MODE ("ok": none).
    task = lsqvce.TaskDefinition(Days=1);
    task.addCommon("c");
    if scheduled
        sys = task.addSystem("S", Arcs=1);
        sys.addLocal("ic", Scope="arc");
    else
        sys = task.addSystem("S");
    end
    sys.addLocal("b", Scope="day");
    sys.addLocal("late", Scope="day", Valid=[5 6]);
    sys.addStream("y", Build=@(day, ctx) tinyBuild(day, ctx, mode, scheduled), Noise=lsqvce.Noise.white(1));
end

function out = tinyBuild(~, ctx, mode, scheduled)
%TINYBUILD  One arc of 30 rows, with the mistake MODE.
    m = 30;
    A = ones(m, 1);
    t = 10 * (0:m-1).';
    l = 2 + 0.1 * cos(t);
    out = ctx.newDay();
    good = {"b", (1:m).' / m};
    if scheduled, good = [{"ic", sin(t)}; good]; end
    switch mode
        case "ok",              out.addArc(A, l, t, Local=good);
        case "unknownName",     out.addArc(A, l, t, Local={"nope", ones(m, 1)});
        case "unknownCols",     ctx.cols("nope");
        case "arcRequired",     ctx.cols("ic");
        case "arcNotScheduled", ctx.cols("ic", Arc=9);
        case "mixed",           ctx.cols(["c" "b"]);
        case "badReturn",       out = 42;
        case "localAndB",       out.addArc(A, l, t, Local=good, B=sparse(m, ctx.NumLocal));
        case "localOneColumn",  out.addArc(A, l, t, Local={"b"});
        case "localWrongSize",  out.addArc(A, l, t, Local={"b", ones(m, 2)});
        case "localTwice",      out.addArc(A, l, t, Local={"b", ones(m, 1); "b", ones(m, 1)});
        case "localNotText",    out.addArc(A, l, t, Local={3, ones(m, 1)});
        case "arcUnscheduled",  out.addArc(A, l, t, Local=good, Arc=4);
        case "arcTwice"
            out.addArc(A, l, t, Local=good, Arc=1);
            out.addArc(A, l, t, Local=good, Arc=1);
        case "arcWithoutSchedule", out.addArc(A, l, t, Local=good, Arc=1);
        case "dayNotOnDay",     out.addArc(A, l, t, Local={"late", ones(m, 1)});
    end
end
