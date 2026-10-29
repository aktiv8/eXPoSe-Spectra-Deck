'use strict';
// Unit tests of the pure half of viewer/viewer.js, run by test_htmlbrowser.py:
//   node tests/viewer_test.js <fixture.json>
// The fixture is written by Python and holds expected values computed there,
// so the JavaScript is checked against the desktop app's own results.
const fs = require('fs');
const path = require('path');
const V = require(path.join(__dirname, '..', 'viewer', 'viewer.js'));
const fx = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));

let failed = 0, ran = 0;
function check(cond, msg) {
  ran++;
  if (!cond) { failed++; console.log('FAIL ' + msg); }
}
function eq(a, b, msg) {
  check(JSON.stringify(a) === JSON.stringify(b), msg + ': ' + JSON.stringify(a) + ' != ' + JSON.stringify(b));
}
function near(a, b, msg, tol) {
  check(Math.abs(a - b) <= (tol || 1e-9), msg + ': ' + a + ' vs ' + b);
}

(async () => {
  // ---- the data Python wrote decodes to the same structure ----
  const data = await V.decode(fx.payload_b64);
  eq(data.v, 1, 'payload version');
  eq(data.samples.length, fx.n_samples, 'sample count');
  const specs = V.prepare(data);
  eq(specs.length, fx.n_spectra, 'spectrum count');
  const s0 = specs[0];
  eq(s0.x.length, fx.first_len, 'axis length');
  near(s0.x[0], fx.first_x0, 'first energy', 1e-6);
  near(s0.x[s0.x.length - 1], fx.first_x_last, 'last energy', 1e-5);
  near(s0.y[3], fx.first_y3, 'a count value', Math.abs(fx.first_y3) * 1e-6 + 1e-12);
  eq(s0.fileName, fx.first_file, 'file name attached');

  // ---- bad input is rejected, not swallowed ----
  let threw = false;
  try { await V.decode('AAAA'); } catch (e) { threw = true; }
  check(threw, 'garbage payload is an error');

  // ---- axes ----
  eq(V.axisValues({ x0: 10, dx: 0.5, n: 4 }), [10, 10.5, 11, 11.5], 'regular axis');
  eq(V.axisValues([1, 2, 4]), [1, 2, 4], 'irregular axis');
  const reg = { binding: true, hv: 1486.6, elabel: 'Binding Energy', eunits: 'eV' };
  let a = V.energyAxis(reg, [285, 290], 'Binding');
  check(a.invert && a.label === 'Binding energy', 'binding axis is inverted');
  a = V.energyAxis(reg, [285, 290], 'Kinetic');
  near(a.x[0], 1201.6, 'KE = hv - BE', 1e-9);
  check(!a.invert && a.label === 'Kinetic energy', 'kinetic axis not inverted');
  a = V.energyAxis({ binding: true, hv: null, elabel: 'Binding Energy' }, [1, 2], 'Kinetic');
  check(!a.ok && a.invert, 'no photon energy: falls back to binding');
  a = V.energyAxis({ binding: false, hv: null, elabel: 'Kinetic Energy' }, [1, 2], 'Kinetic');
  check(a.ok && !a.invert, 'native kinetic axis kept');

  // ---- ion scattering: kinetic by default, energy ratio on request ----
  const idata = await V.decode(fx.iss.payload_b64), ispec = V.prepare(idata)[0];
  near(ispec.reg.iss.e0, fx.iss.e0, 'the payload carries the beam energy', 1e-9);
  a = V.energyAxis(ispec.reg, ispec.x, 'Binding', 'kinetic');
  check(a.ok && !a.invert && !a.ratio && a.label === 'Kinetic energy', 'ISS is kinetic by default');
  eq(a.x, fx.iss.x, 'default ISS axis is the native kinetic one');
  a = V.energyAxis(ispec.reg, ispec.x, 'Binding');
  check(!a.ratio, 'no ISS axis given: kinetic');
  a = V.energyAxis(ispec.reg, ispec.x, 'Binding', 'ratio');
  check(a.ok && a.ratio && !a.invert && a.label === 'Energy ratio' && a.units === 'E/E₀', 'ratio axis look');
  eq(a.x.length, fx.iss.ratio.length, 'ratio axis length');
  a.x.forEach(function (v, i) { near(v, fx.iss.ratio[i], 'ratio value ' + i, 1e-12); });
  near(a.e0, fx.iss.e0, 'the axis keeps the beam energy', 1e-9);
  near(V.markerX({ be: fx.iss.marker_ke, kin: true, hv: null }, a), fx.iss.marker_ratio, 'a kinetic marker at KE / E0', 1e-12);
  eq(V.markerX({ be: 285, kin: false, hv: null }, a), null, 'a binding marker has no place on a ratio axis');
  eq(V.readoutParts(fx.iss.read_v, a, ispec.reg), ['E/E₀ 0.9000', 'KE ' + fx.iss.read_ke.toFixed(1) + ' eV'], 'ratio read-out');
  const noE0 = { binding: false, hv: null, elabel: 'Kinetic Energy', iss: { e0: null } };
  a = V.energyAxis(noE0, [100, 200], 'Binding', 'ratio');
  check(!a.ok && !a.ratio && a.label === 'Kinetic energy', 'no beam energy: kinetic, flagged');
  a = V.energyAxis(reg, [285, 290], 'Binding', 'ratio');
  check(a.ok && a.invert && a.label === 'Binding energy', 'XPS is not affected by the ISS axis');
  a = V.energyAxis(reg, [285, 290], 'Kinetic', 'ratio');
  near(a.x[0], 1201.6, 'XPS kinetic still works beside it', 1e-9);
  // the old markers and read-out are unchanged
  near(V.markerX({ be: 285, kin: false, hv: 1486.6 }, V.energyAxis(reg, [285], 'Binding')), 285, 'binding marker on a binding axis');
  near(V.markerX({ be: 285, kin: false, hv: 1486.6 }, V.energyAxis(reg, [285], 'Kinetic')), 1201.6, 'binding marker on a kinetic axis', 1e-9);
  eq(V.readoutParts(285, V.energyAxis(reg, [285], 'Binding'), reg), ['BE 285.00 eV', 'KE 1201.60 eV'], 'binding read-out');

  // ---- ticks ----
  eq(V.niceTicks(280, 296, 8), [280, 282, 284, 286, 288, 290, 292, 294, 296], 'ticks 280-296');
  eq(V.niceTicks(0, 1, 5), [0, 0.2, 0.4, 0.6, 0.8, 1], 'ticks 0-1');
  eq(V.niceTicks(1, 1, 5), [], 'no ticks for an empty range');
  eq(V.niceTicks(NaN, 1, 5), [], 'no ticks for NaN');
  eq(V.decimals(2), 0, 'decimals of 2');
  eq(V.decimals(0.5), 1, 'decimals of 0.5');
  eq(V.decimals(0.05), 2, 'decimals of 0.05');
  eq(V.niceFloor(1234), 1000, 'niceFloor 1234');
  eq(V.niceFloor(0.37), 0.2, 'niceFloor 0.37');
  eq(V.niceFloor(0), 1, 'niceFloor 0');
  eq(V.fmtY(0), '0', 'fmtY 0');
  eq(V.fmtY(123456), '1.2e5', 'fmtY big');
  eq(V.fmtY(1234.5), '1235', 'fmtY 1234.5');

  // ---- normalising agrees with the desktop app ----
  fx.norm.forEach((c) => {
    near(V.normFactor(c.y, 'max'), c.max, 'max factor ' + JSON.stringify(c.y));
    near(V.normFactor(c.y, 'area'), c.area, 'area factor ' + JSON.stringify(c.y));
    eq(V.normFactor(c.y, 'none'), 1, 'no normalisation');
  });
  fx.range_norm.forEach((c, i) => {
    const a = V.normAffine(c.x, c.y, 'range', c.rng);
    near(a.off, c.off, 'range offset ' + i);
    near(a.div, c.div, 'range divisor ' + i);
  });
  eq(V.normAffine([1, 2], [3, 4], 'none').div, 1, 'range: other modes keep offset 0');

  // ---- colours agree with themes.ramp ----
  fx.ramps.forEach((c) => {
    eq(V.ramp(c.colour, c.n, c.bg), c.expect, 'ramp ' + c.n);
  });
  const cyc = ['#111111', '#222222', '#333333'];
  eq(V.stackColours([0, 0, 0], cyc, '#FFFFFF').length, 3, 'shared slot -> ramp');
  check(V.stackColours([0, 0, 0], cyc, '#FFFFFF')[2] !== '#111111', 'ramp fades');
  eq(V.stackColours([0, 1, 4], cyc, '#FFFFFF'), ['#111111', '#222222', '#222222'], 'categorical wraps');
  eq(V.stackColours([2], cyc, '#FFFFFF'), ['#333333'], 'single trace keeps its colour');

  // ---- labels ----
  const ys = [100, 104, 102, 300];
  const d = V.dodge(ys, 10);
  check(d[0] === 100 && d[2] === 110 && d[1] === 120, 'dodge separates and keeps order: ' + d);
  eq(d[3], 300, 'dodge leaves far labels');
  eq(V.traceLabel({ reg: { level: 3, etch: 90 }, sample: { name: 'S' }, name: 'C 1s', fileName: '' }), 'L3 (90 s)', 'level label');
  eq(V.traceLabel({ reg: { level: null }, sample: { name: 'Sample A' }, name: 'C 1s', fileName: '' }), 'Sample A', 'sample label');
  eq(V.traceLabel({ reg: { level: null }, sample: { name: '' }, name: 'C 1s', fileName: 'run1.vms' }, true), 'run1', 'file label');
  eq(V.traceLabel({ reg: { level: null }, sample: { name: 'x'.repeat(40) }, name: 'C 1s', fileName: '' }).length, 23, 'long label cut');

  // ---- grouping ----
  const g = V.groupSpectra([{ name: 'C 1s' }, { name: 'O 1s' }, { name: 'c  1s' }]);
  eq(g.map((x) => x.items.length), [2, 1], 'grouped by normalised name');
  eq(g.map((x) => x.name), ['C 1s', 'O 1s'], 'first appearance order');

  const og = V.groupSpectra([{ name: 'O 1s' }, { name: 'C 1s' }, { name: 'N 1s' }]);
  eq(V.orderGroups(og, ['c 1s', 'x']).map((x) => x.key), ['c 1s', 'o 1s', 'n 1s'], 'panel order: listed first, rest kept');
  eq(V.moveKey(['a', 'b', 'c'], 'c', 'earlier'), ['a', 'c', 'b'], 'move earlier');
  eq(V.moveKey(['a', 'b', 'c'], 'c', 'onto', 'a'), ['c', 'a', 'b'], 'move onto');
  eq(V.moveKey(['a', 'b', 'c'], 'a', 'last'), ['b', 'c', 'a'], 'move last');

  // ---- interpolation, either direction ----
  near(V.interp([1, 2, 3], [10, 20, 30], 2.5), 25, 'interp ascending');
  near(V.interp([3, 2, 1], [30, 20, 10], 2.5), 25, 'interp descending');
  eq(V.interp([1, 2, 3], [10, 20, 30], 5), null, 'outside the range');
  eq(V.interp([], [], 1), null, 'empty');

  // ---- CSV equals the app's own export ----
  const sel = specs.filter((s) => fx.csv_ids.indexOf(s.id) >= 0);
  const csv = V.buildCsv(sel).replace(/\r\n/g, '\n');
  const want = fx.csv.replace(/\r\n/g, '\n');
  const cl = csv.trim().split('\n'), wl = want.trim().split('\n');
  eq(cl[0], wl[0], 'CSV header equals the app export');
  eq(cl.length, wl.length, 'CSV row count');
  let worst = 0;
  for (let i = 1; i < wl.length; i++) {
    const c = cl[i].split(','), w = wl[i].split(',');
    eq(c.length, w.length, 'CSV columns row ' + i);
    for (let k = 0; k < w.length; k++) {
      if (w[k] === '' || c[k] === '') { eq(c[k], w[k], 'CSV blank cell'); continue; }
      const dv = Math.abs(parseFloat(c[k]) - parseFloat(w[k])) / (Math.abs(parseFloat(w[k])) + 1e-9);
      if (dv > worst) worst = dv;
    }
  }
  check(worst < 1e-6, 'CSV values agree (worst relative error ' + worst + ')');
  eq(V.csvField('a,b'), '"a,b"', 'csv comma quoting');
  eq(V.csvField('say "hi"'), '"say ""hi"""', 'csv quote escaping');
  eq(V.csvField('plain'), 'plain', 'csv plain');
  const lv = V.buildCsv([{ sample: { name: 'S' }, name: 'C 1s', reg: { level: 2, elabel: 'Binding Energy', eunits: 'eV', ylabel: 'Intensity', yunits: 'counts' }, x: [1], y: [2] }]);
  check(lv.split('\r\n')[0].indexOf('S C 1s L2 ') === 0, 'depth levels are told apart in CSV headers');

  // ---- CasaXPS fits: the CSV (with its fit columns) equals the app's export ----
  if (fx.fit) {
    const fdata = await V.decode(fx.fit.payload_b64);
    const fspecs = V.prepare(fdata);
    const fcsv = V.buildCsv(fspecs).replace(/\r\n/g, '\n').trim().split('\n');
    const fwant = fx.fit.csv.replace(/\r\n/g, '\n').trim().split('\n');
    eq(fcsv[0], fwant[0], 'fit CSV header equals the app export (fit columns named the same)');
    eq(fcsv.length, fwant.length, 'fit CSV row count');
    let fworst = 0, blanks = true;
    for (let i = 1; i < fwant.length; i++) {
      const c = fcsv[i].split(','), w = fwant[i].split(',');
      eq(c.length, w.length, 'fit CSV columns row ' + i);
      for (let k = 0; k < w.length; k++) {
        if (w[k] === '' || c[k] === '') { if (c[k] !== w[k]) blanks = false; continue; }
        const dv = Math.abs(parseFloat(c[k]) - parseFloat(w[k])) / (Math.abs(parseFloat(w[k])) + 1e-9);
        if (dv > fworst) fworst = dv;
      }
    }
    check(blanks, 'fit CSV: cells outside the region are blank in both');
    check(fworst < 1e-5, 'fit CSV values agree (worst relative error ' + fworst + ')');
    // curves on every point, null outside the region
    const cur = { i0: 2 }, full = V.curveFull(cur, [5, 6], 6);
    eq(full, [null, null, 5, 6, null, null], 'curveFull places a curve at its first point');
    eq(V.curveFull(null, [1], 3), null, 'no curves, no curve');
    eq(V.residual([10, 20, 30], [9, null, 33]), [1, null, -3], 'residual is data - envelope');
    // colour = the CasaXPS INDEX, computed by casafit.colour_slots on the same cases
    for (const c of fx.fit.slot_cases) {
      eq(V.fitSlots(c.indexes, c.n), c.slots, 'fitSlots(' + JSON.stringify(c.indexes) + ', ' + c.n + ') matches casafit.colour_slots');
    }
    // merging the components of an INDEX group: the same sums as casafit.drawn_components
    for (const c of fx.fit.merge_cases) {
      const comps = c.indexes.map((i) => ({ index: i }));
      const got = V.drawnCurves(comps, c.curves, true, null);
      eq(got.map((e) => ({ j: e.j, members: e.members, values: e.values })), c.expect,
         'drawnCurves merges ' + JSON.stringify(c.indexes) + ' like casafit.drawn_components');
      const plain = V.drawnCurves(comps, c.curves, false, null);
      eq(plain.map((e) => e.values), c.curves, 'without merge every component is drawn on its own');
      eq(plain.map((e) => e.members.length), c.curves.map(() => 1), 'one member each without merge');
    }
    {
      const comps = [{ index: 0 }, { index: 0 }, { index: 0 }, { index: -1 }];
      const curves = [[1, 1], [2, 2], null, [5, 5]];
      eq(V.drawnCurves(comps, curves, true, (j) => j === 1).map((e) => e.values), [[1, 1], [5, 5]],
         'a component the user unticked is left out of the merged sum');
      eq(V.drawnCurves(comps, curves, true, null).map((e) => e.members), [[0, 1], [3]],
         'a component with no curve is skipped, an ungrouped one stands alone');
    }
    // the legend of a stack of fits: one entry per swatch (slot) and name across the traces
    {
      const mk = (comps) => ({ fit: { rows: [{ region: 'C 1s', components: comps }] } });
      const a = mk([{ index: -1, state: 'C 1s c0' }]);                    // ungrouped, slot 0
      const b = mk([{ index: 0, state: 'C 1s c0' }]);                     // INDEX 0, slot 0: the same swatch and words
      const c = mk([{ index: 1, state: 'Oxide' }, { index: -1, state: 'C 1s c0' }]);
      eq(V.fitStatesOf([a, b], 6).map((x) => x.slot + '|' + x.name), ['0|C 1s c0'],
         'the same colour and name in another trace is listed once');
      eq(V.fitStatesOf([a, c], 6).map((x) => x.slot + '|' + x.name), ['0|C 1s c0', '1|Oxide'],
         'a different colour or name is a new entry, in order of appearance');
      eq(V.fitStatesOf([{}, a], 6).length, 1, 'a spectrum without a fit adds nothing');
      eq(V.fitStatesOf([], 6), [], 'no spectra, no states');
      check(V.FIT_STACK_MAX === 8, 'fits on stacks are drawn for up to 8 traces, like the desktop');
      check(V.FIT_LAYERS.indexOf('stack') >= 0 && V.FIT_LAYERS.indexOf('merge') >= 0, 'both fit options are wired');
    }
    for (const n of Object.keys(fx.fit.state_slots)) {
      eq(V.fitStates(fspecs[0].reg, +n).map((x) => x.name + '|' + x.slot),
         fx.fit.state_slots[n].map((x) => x.name + '|' + x.slot), 'chemical states and their INDEX slots for a palette of ' + n);
    }
    eq(V.fitRows({}).length, 0, 'a spectrum without a fit has no rows');
  }

  // ---- figure headings: the page and reportspec.figure_title agree ----
  if (fx.figure_titles) {
    for (const c of fx.figure_titles) {
      eq(V.figureTitle(c.number, c.name, c.sep), c.want, 'figure heading for ' + JSON.stringify(c.name));
    }
    eq(V.figureTitle(2, null), 'Figure 2', 'no name gives the figure number');
    eq(V.figureTitle(2, 'Depth'), 'Figure 2 – Depth', 'the page separates with an en dash by default');
  }

  // ---- the CasaXPS-curves switch picks one of the two row sets everywhere ----
  {
    const rg = { fit: { rows: ['recon'], csv_rows: ['casa'], notes: ['n'], csv_notes: ['cn'] } };
    const plain = { fit: { rows: ['recon'], notes: ['n'] } };
    V.fitSource.csv = false;
    eq(V.fitRows(rg), ['recon'], 'reconstruction rows when the switch is off');
    eq(V.fitNotes(rg), ['n'], 'reconstruction notes when the switch is off');
    V.fitSource.csv = true;
    eq(V.fitRows(rg), ['casa'], 'CasaXPS rows when the switch is on');
    eq(V.fitNotes(rg), ['cn'], 'CasaXPS notes when the switch is on');
    eq(V.fitRows(plain), ['recon'], 'a fit with no CasaXPS match keeps its one set');
    eq(V.fitNotes(plain), ['n'], 'and its notes');
    V.fitSource.csv = false;
  }

  // ---- CasaXPS's own percentages, and a survey that is its own total ----
  {
    const rows = [{ region: 'A', rsf: 1, area: 5, casa_pct: 20 }, { region: 'B', rsf: 1, area: 5, casa_pct: 60 },
                  { region: 'C', rsf: 1, area: 5, casa_why: 'no CasaXPS quantification found for this region' }];
    V.quantSource.casa = true;
    let r = V.quantNormalise(rows, [true, true, true]);
    eq([r[0].at, r[1].at, r[2].at], [25, 75, null], 'CasaXPS numbers are shared out over the counted rows');
    eq(r[2].why, 'no CasaXPS quantification found for this region', 'a region the file lacks is left out and says so');
    r = V.quantNormalise(rows, [true, false, true]);
    eq(r[0].at, 100, 'unticking a region renormalises the rest');
    V.quantSource.casa = false;
    r = V.quantNormalise(rows, [true, true, true]);
    check(Math.abs(r[0].at - 100 / 3) < 1e-9 && r[2].why === '', 'recomputing from the fits ignores the CasaXPS tags');
    V.quantSource.casa = true;
    const reg = { level: null, etch: null, fit: { rows: [{ region: 'C 1s', source: 'survey', components: [] },
                                                          { region: 'C 1s', source: 'high-res', components: [] }] } };
    const g = V.quantGroups([{ id: 's0r0', name: 'S', sample: { id: 's0', name: 'PtCl2' }, reg: reg }]);
    eq(g.map((x) => x.sample + '|' + x.kind), ['PtCl2 (survey)|survey', 'PtCl2|regions'], 'a survey is its own group');
  }

  // ---- fit components: the visible floor, the legend's peak positions ----
  {
    eq(V.clipFloor([0, 5, 100, 2, 0.5, null]), [null, 5, 100, 2, null, null], 'a component stops below 1 % of its own peak');
    eq(V.clipFloor([-100, -0.5, -3]), [-100, null, -3], 'the floor works on the size of a negative component');
    eq(V.clipFloor([null, null]), [null, null], 'a component with no values stays empty');
    eq(V.peakAt([10, 11, 12, 13], [1, 9, 9, null]), 11, 'the first of equal maxima');
    eq(V.peakAt([1, 2], [null, null]), null, 'no peak without values');
    const st = [{ key: 'a|i0', name: 'Pt 4f' }, { key: 'a|i1', name: 'Pt 4f' }, { key: 'a|i2', name: 'Cl 2p' },
                { key: 'a|c3', name: 'A very long state name indeed' }];
    eq(V.legendLabels(st, { 'a|i0': 71.04, 'a|i1': 74.4 }),
       ['Pt 4f (71.0)', 'Pt 4f (74.4)', 'Cl 2p', 'A very long state nam\u2026'], 'states that share a name show where they peak');
    eq(V.legendLabels(st.slice(0, 2), {}), ['Pt 4f', 'Pt 4f'], 'no peaks known: the plain name');
    eq(V.stateKey({ region: 'Pt 4F', components: [{ index: 2 }, { index: -1 }] }, 1), 'pt 4f|c1', 'an ungrouped state is keyed by position');
  }

  // ---- camera pictures: zoom and pan stay inside the picture ----
  {
    const W = 1600, H = 800;
    let v = V.camZoom(null, W, H, 0.5, 800, 400);
    eq(v, { x: [400, 1200], y: [200, 600] }, 'zooming in about the centre halves the window');
    v = V.camZoom(null, W, H, 0.5, 0, 0);
    eq(v, { x: [0, 800], y: [0, 400] }, 'zooming about a corner keeps that corner');
    eq(V.camZoom(v, W, H, 2, 0, 0), null, 'zooming back out gives the whole picture');
    eq(V.camPan(v, W, H, 5000, 5000), { x: [800, 1600], y: [400, 800] }, 'a pan stops at the edge of the picture');
    eq(V.camPan(null, W, H, 10, 10), null, 'there is nothing to pan on the whole picture');
    let z = null;
    for (let i = 0; i < 40; i++) z = V.camZoom(z, W, H, 0.5, 800, 400);
    check(z.x[1] - z.x[0] >= W / V.CAM_MAX_ZOOM - 1e-9, 'the zoom has a limit');
    eq(V.camFit([-50, 750], [-20, 380], W, H), { x: [0, 800], y: [0, 400] }, 'a window past the corner slides back in');
  }

  // ---- zoom stays inside the data ----
  {
    eq(V.clampView(525, 535, 518, 547), [525, 535], 'a zoom inside the data is kept');
    eq(V.clampView(505, 520, 518, 547), [518, 533], 'a zoom past the low edge slides back in');
    eq(V.clampView(535, 555, 518, 547), [527, 547], 'a zoom past the high edge slides back in');
    eq(V.clampView(0, 2000, 518, 547), null, 'a zoom wider than the data is the full view');
    eq(V.clampView(518, 547, 518, 547), null, 'a zoom as wide as the data is the full view');
    eq(V.clampView(1000, 1010, 518, 547), null, 'a zoom wholly outside the data is dropped');
    eq(V.clampView(530, 530, 518, 547), null, 'an empty zoom is dropped');
    eq(V.clampView(NaN, 530, 518, 547), null, 'a zoom that is not a number is dropped');
  }

  // ---- quantification: the same table as quant.csv_rows ----
  if (fx.fit && fx.fit.quant) {
    const qx = fx.fit.quant;
    const cmp = (got, want, label) => {
      eq(got.length, want.length, label + ' row count');
      got.forEach((row, i) => {
        eq(row.length, want[i].length, label + ' columns row ' + i);
        row.forEach((cell, k) => {
          const w = want[i][k], numeric = i > 0 && k >= 5 && k <= 8 || i > 0 && k === 10;
          if (numeric && cell !== '' && w !== '') near(parseFloat(cell), parseFloat(w), label + ' cell ' + i + ',' + k, Math.abs(parseFloat(w)) * 1e-5 + 1e-9);
          else eq(cell, w, label + ' cell ' + i + ',' + k);
        });
      });
    };
    cmp(V.quantTable(qx.groups, null, false), qx.plain, 'quant table');
    cmp(V.quantTable(qx.groups, qx.excluded_keys, false), qx.excluded, 'quant table with a region left out');
    cmp(V.quantTable(qx.groups, null, true), qx.trans, 'quant table with the transmission function');
    const csv = V.quantCsv(qx.groups, null, false).split('\r\n');
    eq(csv[0], V.QUANT_HEADER.map(V.csvField).join(','), 'quant CSV header');
    check(csv[1].indexOf('"S, 1",1,C 1s,C 1s,Shirley,0.278') === 0, 'quant CSV quotes a comma in the sample name: ' + csv[1]);
    const n = V.quantNormalise([{ area: 1, rsf: 1 }, { area: 3, rsf: 1 }], null, false);
    near(n[0].at, 25, 'atomic percent'); near(n[1].at, 75, 'atomic percent');
    eq(V.quantNormalise([{ area: 5, rsf: 0 }], null, false)[0].why, 'no RSF', 'no RSF is said');

    // ---- RSF fallback tiers (mirrors quant.py's own tiers, test_quant.py) ----
    {
      // tier: a component's own RSF, summed, when the region has none --
      // Al2O3.kfit's real Al 2p shape (region rsf 0, components 0.37/0.19)
      const r = V.quantNormalise([{ area: 999, rsf: 0, components: [
        { area: 111265.21, rsf: 0.37 }, { area: 54535.04, rsf: 0.19 }] }], null, false)[0];
      near(r.corrected, 111265.21 / 0.37 + 54535.04 / 0.19, 'component-tier corrected');
      eq(r.rsfSource, 'component', 'component-tier rsfSource');

      // tier: no usable component RSF falls through to "no RSF"
      eq(V.quantNormalise([{ area: 10, rsf: 0, components: [{ area: 5, rsf: 0 }] }],
        null, false)[0].why, 'no RSF', 'no usable component RSF is still no RSF');

      // tier: the reference table, off by default, on when supplied
      const table = [{ 0: 'scofield', 1: 'Al', 2: 'Pt 4f', 3: 15.45 }].map((o) => [o[0], o[1], o[2], o[3]]);
      const ptRow = { area: 100, rsf: 0, region: 'Pt 4f', photon_energy: 1486.6, components: [] };
      eq(V.quantNormalise([ptRow], null, false)[0].why, 'no RSF', 'table fallback is off by default');
      const t = V.quantNormalise([ptRow], null, false, table, 'scofield')[0];
      near(t.corrected, 100 / 15.45, 'table-tier corrected');
      eq([t.rsfSource, t.rsfAnode, t.rsfValue], ['scofield', 'Al', 15.45], 'table-tier provenance');

      // a component-level fix always wins over the table (the file's own data)
      const compRow = { area: 10, rsf: 0, region: 'Pt 4f', photon_energy: 1486.6,
                       components: [{ area: 5, rsf: 2 }] };
      eq(V.quantNormalise([compRow], null, false, table, 'scofield')[0].rsfSource, 'component',
        'component tier beats the table fallback');

      // tier: scofield_tpp2m/scofield_ke06 (mirrors test_quant.py's own cases)
      const mfpRow = { area: 100, rsf: 0, region: 'Pt 4f', photon_energy: 1486.6,
                      be_lo: 70, be_hi: 80, components: [] };
      const ke = 1486.6 - 75;
      const tpp2m = V.quantNormalise([mfpRow], null, false, table, 'scofield_tpp2m')[0];
      eq(tpp2m.rsfSource, 'scofield_tpp2m', 'scofield_tpp2m rsfSource');
      near(tpp2m.imfpNm, V.imfpNm(ke), 'scofield_tpp2m imfpNm matches V.imfpNm');
      near(tpp2m.rsfValue, 15.45 * V.imfpNm(ke), 'scofield_tpp2m rsfValue');
      eq(tpp2m.kePowerFactor, null, 'scofield_tpp2m leaves kePowerFactor null');

      // CasaXPS's own TPP-2M matrix: same value as imfp.casa_imfp_nm (Python)
      // and no energy limit
      const casa = V.quantNormalise([mfpRow], null, false, table, 'scofield_tpp2m_casa')[0];
      eq(casa.rsfSource, 'scofield_tpp2m_casa', 'scofield_tpp2m_casa rsfSource');
      near(casa.imfpNm, V.imfpNmCasa(ke), 'scofield_tpp2m_casa imfpNm matches V.imfpNmCasa');
      near(casa.imfpNm, 3.6732544441977097, 'Casa TPP-2M equals imfp.casa_imfp_nm(1411.6)', 1e-9);
      const casaLow = V.quantNormalise([{ area: 100, rsf: 0, region: 'Pt 4f', photon_energy: 1486.6,
                      be_lo: 1440, be_hi: 1442, components: [] }], null, false, table, 'scofield_tpp2m_casa')[0];
      eq(casaLow.rsfSource, 'scofield_tpp2m_casa', 'Casa TPP-2M has no energy floor');

      const ke06 = V.quantNormalise([mfpRow], null, false, table, 'scofield_ke06')[0];
      eq(ke06.rsfSource, 'scofield_ke06', 'scofield_ke06 rsfSource');
      near(ke06.kePowerFactor, V.kePowerFactor(ke), 'scofield_ke06 kePowerFactor matches V.kePowerFactor');
      near(ke06.rsfValue, 15.45 * V.kePowerFactor(ke), 'scofield_ke06 rsfValue');
      eq(ke06.imfpNm, null, 'scofield_ke06 leaves imfpNm null');

      // out of TPP-2M's valid range (or missing be_lo/be_hi): falls back to
      // the plain, unmultiplied Scofield value, not a silent extrapolation
      const lowKeRow = { area: 100, rsf: 0, region: 'Pt 4f', photon_energy: 1486.6,
                        be_lo: 1440, be_hi: 1442, components: [] };
      const fallback = V.quantNormalise([lowKeRow], null, false, table, 'scofield_tpp2m')[0];
      eq(fallback.rsfSource, 'scofield', 'out-of-range TPP-2M falls back to plain scofield');
      eq(fallback.imfpNm, null, 'out-of-range TPP-2M leaves imfpNm null');
      near(fallback.rsfValue, 15.45, 'out-of-range TPP-2M value is unmultiplied');

      // the empirical Kratos Axis F1s library must never get an IMFP factor
      const kratosRow = { area: 100, rsf: 0, region: 'Pt 4f', photon_energy: 1486.6,
                         be_lo: 70, be_hi: 80, components: [] };
      const kratosTable = [['kratos_f1s', 'Al', 'Pt 4f', 5.58]];
      const kr = V.quantNormalise([kratosRow], null, false, kratosTable, 'kratos_f1s')[0];
      eq(kr.rsfSource, 'kratos_f1s', 'kratos_f1s rsfSource');
      eq(kr.imfpNm, null, 'kratos_f1s is never IMFP-corrected');
    }

    // ---- imfpNm / kePowerFactor (mirrors imfp.py, tests/test_imfp.py) ----
    {
      // hand-computed from the same formula independently (test_imfp.py's own reference table)
      const REF = { 50: 0.5066507922262123, 100: 0.577623423099999, 200: 0.7947830339179165,
                    500: 1.4263303268788368, 1000: 2.368001423854056, 2000: 4.060332098481516 };
      Object.keys(REF).forEach((k) => near(V.imfpNm(+k), REF[k], 'imfpNm reference KE=' + k, 1e-9));
      eq(V.imfpNm(49.9), null, 'imfpNm below the valid range');
      eq(V.imfpNm(2000.1), null, 'imfpNm above the valid range');
      eq(V.imfpNm(0), null, 'imfpNm with no kinetic energy');
      near(V.kePowerFactor(1419), Math.pow(1419, 0.6), 'kePowerFactor matches KE^0.6');
      eq(V.kePowerFactor(0), null, 'kePowerFactor with no kinetic energy');
      eq(V.kePowerFactor(-5), null, 'kePowerFactor with a negative kinetic energy');
    }

    // ---- anodeFor / rsfOf (mirrors rsf.py) ----
    eq(V.anodeFor(null), 'Al', 'anodeFor defaults to Al');
    eq(V.anodeFor(1486.6), 'Al', 'anodeFor Al Kalpha');
    eq(V.anodeFor(1253.6), 'Mg', 'anodeFor Mg Kalpha');
    eq(V.anodeFor(1486.68), 'Al', 'anodeFor Thermo Al Kalpha');
    eq(V.anodeFor(5414.8), null, 'anodeFor Cr 5414.8 eV (HAXPES) has no anode');
    eq(V.anodeFor(2984.2), null, 'anodeFor Ag Lalpha has no anode');
    {
      const table = [['scofield', 'Al', 'Pt 4f', 15.45], ['scofield', 'Mg', 'Pt 4f', 15.86]];
      eq(V.rsfOf('Pt 4f', table, 'scofield', 5414.8), null, 'rsfOf Cr: no table value');
      near(V.rsfOf('Pt 4f', table, 'scofield', 1486.6), 15.45, 'rsfOf Al');
      near(V.rsfOf('Pt 4f', table, 'scofield', 1253.6), 15.86, 'rsfOf Mg');
      eq(V.rsfOf('Zz 9z', table, 'scofield', 1486.6), null, 'rsfOf no match');
    }

    // ---- elementOf / preferredDefaults (mirrors resultspages.py) ----
    eq(V.elementOf('Cl 2p'), 'Cl', 'elementOf a real region');
    eq(V.elementOf('WideScan'), '', 'elementOf a non-region string');
    {
      // the user's own real case: Pt 4f is preferred over Pt 4d
      const entries = [{ key: 'a', row: { region: 'Pt 4f' } }, { key: 'b', row: { region: 'Pt 4d' } }];
      eq(V.preferredDefaults(entries, {}), { b: false }, 'Pt 4d defaults to unticked, Pt 4f untouched');
      // Ti has no preference entry: nothing is seeded
      const ti = [{ key: 'a', row: { region: 'Ti 2p' } }, { key: 'b', row: { region: 'Ti 1s' } }];
      eq(V.preferredDefaults(ti, {}), {}, 'an element with no preference entry is untouched');
      // an already-excluded entry is not treated as a surviving competitor
      eq(V.preferredDefaults(entries, { a: false }), {}, 'nothing left to prefer when Pt 4f is already unticked');
    }
    {
      // the app's own region ticks (payload.quant_include) start the tick map
      const inc = {};
      V.seedInclude(inc, { 's0r0:0': false, 's0r0:1': true, bad: 'no', worse: 1 });
      eq(inc, { 's0r0:0': false, 's0r0:1': true }, 'only booleans are seeded');
      V.seedInclude(inc, undefined);
      V.seedInclude(inc, null);
      eq(inc, { 's0r0:0': false, 's0r0:1': true }, 'a missing payload key changes nothing');
      // renderFitQuant fills preferredDefaults only where a key is absent, so
      // the app's explicit tick on a non-preferred line wins
      const entries = [{ key: 'a', row: { region: 'Pt 4f' } }, { key: 'b', row: { region: 'Pt 4d' } }];
      const ticks = V.seedInclude({}, { b: true });
      const seed = V.preferredDefaults(entries, ticks);
      Object.keys(seed).forEach((k) => { if (!(k in ticks)) ticks[k] = seed[k]; });
      eq(ticks, { b: true }, 'the app tick on Pt 4d beats the default');
      const unticked = V.seedInclude({}, { a: false });
      const seed2 = V.preferredDefaults(entries, unticked);
      Object.keys(seed2).forEach((k) => { if (!(k in unticked)) unticked[k] = seed2[k]; });
      eq(unticked, { a: false }, 'with Pt 4f unticked nothing else is defaulted off');
    }
    {
      // any element counts from one line: the most sensitive (highest RSF) unless a standard line is fitted
      const r = (region, rsf, area) => ({ region: region, rsf: rsf, area: area, components: [] });
      const ti = [{ key: 'a', row: r('Ti 3p', 0.8, 40) }, { key: 'b', row: r('Ti 2p', 5.22, 100) }, { key: 'c', row: r('O 1s', 2.93, 50) }];
      eq(V.preferredDefaults(ti, {}), { a: false }, 'Ti 2p (higher RSF) wins over Ti 3p');
      eq(V.preferredDefaults(ti, { b: false }), {}, 'nothing to choose between once one line is unticked');
      const none = [{ key: 'a', row: r('Ti 3p', 0, 40) }, { key: 'b', row: r('Ti 2p', 5.22, 100) }];
      eq(V.preferredDefaults(none, {}), {}, 'a line without an RSF is no competitor');
      const pt = [{ key: 'a', row: r('Pt 4d', 19.8, 300) }, { key: 'b', row: r('Pt 4f', 15.45, 200) }];
      eq(V.preferredDefaults(pt, {}), { a: false }, 'the standard line still beats a higher RSF');
      eq(V.siblings(ti, 'a'), ['b'], 'siblings: the same element, another line');
      eq(V.siblings(ti, 'c'), [], 'siblings: none for a lone element');
      const dup = [{ key: 'a', row: r('C 1s', 1, 1) }, { key: 'b', row: r('C 1s', 1, 1) }];
      eq(V.siblings(dup, 'a'), [], 'a second copy of one region is not a sibling');
    }
    {
      // Avantage's Area (N) = Area (P) / (SF x KE^0.6), shared out by state (mirrors quant.comp_term / states)
      const c = (n, p, ke) => ({ gk: n, state: n, area: p * 1.05, area_p: p, rsf: 1.0, ke: ke, ecsf: 0.6 });
      near(V.compTerm(c('A', 10954.67, 1201.745)), 155.49, 'compTerm: Area (N) of the PET C 1s phenyl peak', 0.02);
      near(V.compTerm({ area: 6, rsf: 2, ke: 1000 }), 3, 'compTerm: no exponent is the plain ratio');
      eq(V.compTerm({ area: 6, rsf: 0 }), null, 'compTerm: no RSF');
      const row = { rsf: 0, components: [c('A', 10954.67, 1201.745), c('B', 4379.85, 1200.199), c('C', 3578.5, 1197.795), c('D', 943.5, 1195.459)] };
      eq(V.quantStates(row, 100).map((s) => +s.at.toFixed(2)), [55.13, 22.06, 18.05, 4.76], 'states follow Area (N)');
      eq(V.quantNormalise([row], null, false, null, '')[0].rsfSource, 'component', 'component RSF tier');
      const tp = { area: 1000, area_p: 900, rsf: 2.0, ke: 1200, mfp: 'tpp2m', ecsf: 0 };
      near(V.compTerm(tp), 900 / (2.0 * Math.pow(1118, 0.6) * V.imfpNm(1200) / V.imfpNm(1118)), 'compTerm: TPP-2M term (ALTHERMO1)', 1e-9);
      near(V.compTerm(Object.assign({}, tp, { ke: 3000 })), 450, 'compTerm: outside TPP-2M range, the plain ratio');
      near(V.compTerm(Object.assign({}, tp, { tx: 116.5 })), V.compTerm(tp) / 116.5, 'compTerm: TXFN divided out as well', 1e-9);
      const own = { rsf: 0.3, components: [{ gk: 'a', state: 'A', area: 30, rsf: 0.3, ke: 1000, ecsf: 0.6 }, { gk: 'b', state: 'B', area: 10, rsf: 0.3, ke: 800, ecsf: 0.6 }] };
      near(V.quantStates(own, 100)[0].at, 75, 'a region with its own RSF shares out by area');
    }
    eq(V.quantStates({ components: [{ gk: 'i1', state: 'A', area: 3 }, { gk: 'i1', state: 'A', area: 1 }, { gk: 'nB', state: 'B', area: -2 }] }, 40)
      .map((s) => s.name + ':' + s.at), ['A:40', 'B:0'], 'states share the row, negative areas count as zero');
  }

  // ---- depth profiles: the same series as quant.profile ----
  if (fx.fit && fx.fit.quant && fx.fit.quant.profile) {
    const px = fx.fit.quant.profile;
    const same = (got, want, label) => {
      eq(got.levels, want.levels, label + ' levels');
      eq(got.series.map((s) => s.name), want.series.map((s) => s.name), label + ' series names');
      got.series.forEach((s, k) => s.values.forEach((v, i) => {
        const w = want.series[k].values[i];
        if (v === null || w === null) eq(v, w, label + ' gap ' + s.name + ' ' + i);
        else near(v, w, label + ' ' + s.name + ' ' + i, Math.abs(w) * 1e-9 + 1e-9);
      }));
    };
    ['element', 'state', 'share'].forEach((m) => same(V.profile(px.groups, m, null, false), px[m], 'profile ' + m));
    same(V.profile(px.groups, 'element', px.exclude, false), px.element_excl, 'profile with a region left out');
  }
  {
    // where the levels sit, and which axes they support
    const info = V.levelInfo({ level: 2, etch: 60, meta: { 'Depth (nm)': '3.5', 'Fluence (ions/cm²)': '1.2e+16' } });
    eq([info.level, info.etch, info.depth, info.fluence], [2, 60, 3.5, 1.2e16], 'levelInfo reads the metadata');
    eq(V.levelInfo({ level: 1, meta: {} }).depth, null, 'no depth without sputter settings');
    const lv = (etch, depth) => ({ level: 0, etch: etch, depth: depth, fluence: null });
    eq(V.profileAxes([lv(0, 0), lv(30, 1.5), lv(60, 3)].map((x, i) => Object.assign(x, { level: i }))).map((a) => a.id),
       ['depth', 'etch', 'level'], 'axes with depth, etch time and level');
    eq(V.profileAxes([lv(0, null), lv(0, null), lv(0, null)].map((x, i) => Object.assign(x, { level: i }))).map((a) => a.id),
       ['level'], 'all-zero etch times are not an axis');
    const specs = [
      { name: 'C 1s', reg: { level: 1 }, y: [1, 9, 4] }, { name: 'C 1s', reg: { level: 0 }, y: [2, 3] },
      { name: 'O 1s', reg: { level: 0 }, y: [7] } ];
    const pm = V.peakMaxSeries(specs);
    eq(pm.levels, [0, 1], 'peak maximum levels sorted');
    eq(pm.series, [{ name: 'C 1s', values: [3, 9] }, { name: 'O 1s', values: [7, null] }], 'peak maxima per level');
    const tb = V.profileTable({ levels: [0, 1], series: [{ name: 'A', values: [12.3456789, null] }] },
      [{ level: 0, etch: 0, depth: null, fluence: null }, { level: 1, etch: 30, depth: null, fluence: null }]);
    eq(tb, [['Level', 'Etch time (s)', 'A'], ['0', '0', '12.3457'], ['1', '30', '']], 'profile table');
  }

  // ---- element identification: the same candidates, in the same order, as xpslines.py ----
  if (fx.elements) {
    const ex = fx.elements;
    ex.cases.forEach((c, i) => {
      const got = V.candidates(c.be, c.win, ex.table, c.hv, c.split);
      eq(got.map((x) => x.label), c.labels, 'candidates ' + i + ' (' + c.be + ' eV)');
      got.forEach((x, k) => near(x.d, c.deltas[k], 'candidate delta ' + i + '.' + k, 1e-9));
    });
    check(ex.cases.some((c) => c.labels.length > 1), 'the fixture has an ambiguous peak, so the order is tested');
    check(ex.table.lines.some((l) => l[2] === null), 'the fixture has Auger lines');
    near(V.lineBe(['O', 'KLL', null, 510, 1], 1486.6, 1486.6), 976.6, 'Auger binding energy follows hv', 1e-9);
    near(V.lineBe(['O', 'KLL', null, 510, 1], null, 1253.6), 743.6, 'default hv when none is given', 1e-9);
    eq(V.candidates(700, 1, ex.table, 1486.6, false), [], 'nothing near 700 eV within 1 eV');
    eq(V.baseLabel('Ti 2p3/2'), 'Ti 2p', 'base label of a component');
    eq(V.baseLabel('Ti 2p'), 'Ti 2p', 'base label of a pair');
    eq(V.splitLine('4f7/2'), ['4f', '7/2'], 'split line');
    eq(V.splitLine('KL1'), ['KL1', ''], 'an Auger line does not split');
    eq(V.lineLabel(['C', '1s', 285, null, 1]), 'C 1s', 'line label');
    (ex.nearby || []).forEach((c, i) => {
      const got = V.nearbyLines(c.be, c.win, ex.table, c.hv, c.exclude, c.secondary, c.auger, 2, c.split);
      eq(got.map((x) => x.label), c.labels, 'nearby ' + i + ' labels');
      eq(got.map((x) => x.tier), c.tiers, 'nearby ' + i + ' tiers');
      got.forEach((x, k) => near(x.be, c.candidate_be[k], 'nearby ' + i + '.' + k + ' be', 1e-9));
    });
  }

  // ---- CasaXPS's own exported quantification (casaquant.py) ----
  if (fx.casaxps) {
    const cx = fx.casaxps;
    eq(V.casaxpsSurveyRows(cx.sample_casaxps), cx.survey_rows, 'casaxps survey rows');
    eq(V.casaxpsRegionsRows(cx.sample_casaxps), cx.regions_rows, 'casaxps regions rows');
    eq(V.casaxpsDparamRows(cx.sample_casaxps), cx.dparam_rows, 'casaxps dparam rows');
    const scsv = V.casaxpsCsv(cx.survey_rows).split('\r\n');
    eq(scsv[0], 'Element,%Conc', 'casaxps CSV header');
    eq(scsv[1], 'O 1s,1.82', 'casaxps CSV row');
    const cdata = await V.decode(cx.payload_b64);
    eq(cdata.samples[0].casaxps, cx.sample_casaxps,
       'casaxps travels through the payload encode/decode round trip');
    const cspecs = V.prepare(cdata);
    check(V.quantGroups(cspecs).length === 0, 'the fixture region has no fit: no fit-derived group');
    const bundle = V.bundleFiles(cdata, cspecs);
    const names = bundle.map((f) => f.name);
    check(names.indexOf('casaxps/A_survey.csv') >= 0, 'bundle includes the casaxps survey CSV');
    check(names.indexOf('casaxps/A_regions.csv') >= 0, 'bundle includes the casaxps regions CSV');
    check(names.indexOf('casaxps/A_dparam.csv') >= 0, 'bundle includes the casaxps D-parameter CSV');
    check(bundle[0].data.indexOf('casaxps/') >= 0, 'README mentions the casaxps folder');
  }

  // ---- the ZIP the page hands over: read back by Python's zipfile ----
  check(V.crc32(V.utf8('123456789')) === 0xCBF43926, 'crc32 of the standard check string');
  eq(V.safeName('a/b:c'), 'a_b_c', 'safeName replaces separators');
  eq(V.safeName('..'), 'unnamed', 'safeName never gives a dot name');
  eq(V.safeName('Pt #001a (1)'), 'Pt _001a (1)', 'safeName keeps letters, digits, spaces and brackets');
  const bundle = V.bundleFiles(data, specs);
  eq(bundle.map((f) => f.name), ['README.txt', 'csv/A/C 1s.csv', 'csv/A/O 1s.csv', 'csv/B/C 1s.csv'], 'bundle contents');
  check(bundle[0].data.indexOf('3 spectra from 2 samples') >= 0, 'README counts the spectra');
  check(bundle[1].data.charAt(0) === '﻿', 'CSV files start with a byte-order mark for Excel');
  if (fx.zip_out) {
    fs.writeFileSync(fx.zip_out, Buffer.from(V.zip(bundle.concat([{ name: 'ünï/tëst.txt', data: 'héllo' }]), new Date(2026, 0, 2, 3, 4, 6))));
  }

  // ---- metadata summary ----
  const sm = V.summariseMeta([{ A: '1', B: 'x', C: '' }, { A: '1', B: 'y', C: '' }, { A: '1', B: 'x', D: '5' }], []);
  eq(sm.common, [['A', '1']], 'common metadata');
  eq(sm.varying, ['B', 'D'], 'varying metadata');
  eq(V.summariseMeta([{ Sample: 'a', Q: '1' }, { Sample: 'b', Q: '1' }], ['Sample']).common, [['Q', '1']], 'skipped keys');

  // ---- paragraphs ----
  eq(V.paragraphs('a\nb\n\n\nc\r\n\r\nd\n'), ['a\nb', 'c', 'd'], 'paragraphs');
  eq(V.paragraphs(''), [], 'no paragraphs');

  // ---- SnapMaps: the page's maths against snapmap.py on the same map ----
  if (fx.map) {
    const m = fx.map.entry, e = fx.map.energy;
    const data = await V.decodeMap(m);
    eq(data.length, m.nx * m.ny * m.n, 'map data length');
    eq(V.axisValues(m.e).length, m.n, 'map energy axis length');
    const win = fx.map.win;
    const close = (a, b, msg, tol) => {
      let worst = 0;
      eq(a.length, b.length, msg + ' length');
      for (let i = 0; i < b.length; i++) worst = Math.max(worst, Math.abs(a[i] - b[i]));
      check(worst <= (tol || 1e-6), msg + ' (worst difference ' + worst + ')');
    };
    eq(V.mapChannels(e, fx.map.channels[0][0], fx.map.channels[0][1]), fx.map.channels[1], 'channels in a window');
    eq(V.mapChannels(e, fx.map.channels[0][1], fx.map.channels[0][0]), fx.map.channels[1], 'window given either way round');
    close(V.mapImage(m, data, e, win[0], win[1], false), fx.map.image, 'map image');
    close(V.mapImage(m, data, e, win[0], win[1], true), fx.map.image_bg, 'map image, background removed');
    eq(Array.from(V.mapImage(m, data, e, 1e9, 2e9, false)).every((v) => v === 0), true, 'empty window is zeros');
    close(V.meanSpectrum(m, data, null), fx.map.total_mean, 'whole-map mean spectrum');
    const r = V.rectMask(m, fx.map.rect[0], fx.map.rect[1], fx.map.rect[2], fx.map.rect[3]);
    eq(Array.from(r.mask), fx.map.mask, 'rectangle mask');
    eq(r.count, fx.map.mask_count, 'rectangle pixel count');
    eq(Array.from(V.rectMask(m, fx.map.rect[2], fx.map.rect[3], fx.map.rect[0], fx.map.rect[1]).mask), fx.map.mask, 'rectangle corners either way');
    close(V.meanSpectrum(m, data, r.mask), fx.map.roi_mean, 'area mean spectrum');
    fx.map.pixels.forEach((p) => eq(V.pixelAt(m, p[0], p[1]), p[2], 'pixel at ' + p[0] + ',' + p[1]));
    const rng = V.colourRange(V.mapImage(m, data, e, win[0], win[1], false));
    near(rng[0], fx.map.range[0], 'colour range low', 1e-6);
    near(rng[1], fx.map.range[1], 'colour range high', 1e-6);
    eq(V.colourRange([NaN, Infinity]), [0, 1], 'colour range of nothing');
    eq(V.colourRange([5, 5, 5]), [5, 6], 'colour range of a flat map');
    fx.map.window_cases.forEach((c, i) => {
      const w = V.defaultWindow(c.energy, c.y);
      near(w[0], c.expect[0], 'default window low, case ' + i, 1e-9);
      near(w[1], c.expect[1], 'default window high, case ' + i, 1e-9);
    });
    const csv = V.mapCsv(m, fx.map.image).trim().split('\r\n'), want = fx.map.csv.trim().split('\n');
    eq(csv.length, want.length, 'map CSV rows');
    eq(csv[0], want[0], 'map CSV header');
    let cw = 0;
    for (let i = 1; i < want.length; i++) {
      const a = csv[i].split(','), b = want[i].split(',');
      eq(a.length, b.length, 'map CSV columns row ' + i);
      for (let k = 0; k < b.length; k++) cw = Math.max(cw, Math.abs(parseFloat(a[k]) - parseFloat(b[k])) / (Math.abs(parseFloat(b[k])) + 1e-9));
    }
    check(cw < 1e-5, 'map CSV values agree (worst relative ' + cw + ')');
    eq(V.scaleColour('Viridis', 0).map(Math.round), [68, 1, 84], 'viridis starts dark purple');
    eq(V.scaleColour('Viridis', 1).map(Math.round), [253, 231, 37], 'viridis ends yellow');
    eq(V.scaleColour('Nope', 0.5).length, 3, 'unknown scale falls back');
    eq(V.scaleColour('Greys', -5).map(Math.round), V.scaleColour('Greys', 0).map(Math.round), 'colour position clamps');
    let bad = false;
    try { await V.decodeMap(Object.assign({}, m, { n: m.n + 1 })); } catch (err) { bad = true; }
    check(bad, 'a map of the wrong size is an error');
  }

  // ---- Kratos imaging maps: the page's maths against kratosmap.py on the same maps ----
  if (fx.imaging) {
    const im = fx.imaging, d2 = await V.decode(im.payload_b64), maps = d2.imaging;
    eq(maps.length, im.counts.length, 'image map count');
    eq(maps[0].n, 1, 'one channel');
    const counts = [];
    for (const m of maps) counts.push(V.imageCounts(m, await V.decodeMap(m)));
    counts.forEach((c, i) => eq(Array.from(c), im.counts[i], 'counts of map ' + i + ' are exactly the recorded ones'));
    const m0 = maps[0], nx = m0.nx, ny = m0.ny;
    Object.keys(im.blur).forEach((sg) => {
      const b = V.imageBlur(counts[0], nx, ny, +sg);
      let worst = 0;
      for (let i = 0; i < b.length; i++) worst = Math.max(worst, Math.abs(b[i] - im.blur[sg][i]));
      check(worst < 1e-9, 'blur sigma ' + sg + ' (worst ' + worst + ')');
    });
    Object.keys(im.focus).forEach((sg) => counts.forEach((c, i) => {
      near(V.imageFocus(c, nx, ny, +sg), im.focus[sg][i], 'sharpness of map ' + i + ' at sigma ' + sg, 1e-9);
    }));
    eq(V.imageFocus(new Float64Array(nx * ny), nx, ny), 0, 'sharpness of an empty image');
    const arr = (a) => a.map((v) => (v === null ? NaN : v));
    const mean = V.imageRoiMeans(counts, Uint8Array.from(im.mask));
    mean.forEach((v, i) => near(v, im.roi_means[i], 'area mean of map ' + i, 1e-9));
    V.imageRoiMeans(counts, null).forEach((v, i) => near(v, im.whole_means[i], 'whole-image mean of map ' + i, 1e-9));
    check(Number.isNaN(V.imageRoiMeans(counts, new Uint8Array(3))[0]), 'a mask of another size gives NaN');
    check(Number.isNaN(V.imageRoiMeans(counts, new Uint8Array(nx * ny))[0]), 'an empty mask gives NaN');
    Object.keys(im.select).forEach((k) => {
      const [cur, how] = k.split('|');
      eq(V.imageSelect(maps, maps[+cur], how).map((m) => maps.indexOf(m)), im.select[k], 'frames shown for ' + k);
    });
    maps.forEach((m, i) => eq(V.imageDefaultFilter(maps, m), im.default[i], 'default filter of map ' + i));
    Object.keys(im.axes).forEach((k) => {
      const a = im.axes[k], got = V.imageSeriesAxis(a.idx.map((i) => maps[i]));
      eq(got[0], a.label, 'series axis label ' + k);
      got[1].forEach((v, j) => near(v, a.x[j], 'series axis value ' + k + ' ' + j, 1e-9));
    });
    const csv = V.mapCsv(m0, counts[0]).trim().split('\r\n'), want = im.csv.trim().split('\n');
    eq(csv.length, want.length, 'image CSV rows');
    eq(csv[0], want[0], 'image CSV header');
    for (let i = 1; i < want.length; i++) eq(csv[i].split(',').map(Number), want[i].split(',').map(Number), 'image CSV row ' + i);
    const table = V.imageTableCsv(maps, mean, mean.map((_v, i) => i)).trim().split('\r\n');
    eq(table.length, maps.length + 1, 'table CSV rows');
    eq(table[0].split(',').length, 10, 'table CSV columns');
    eq(table[1].split(',').slice(0, 4), ['1', 'Au 4f', 'P1', '83.91'], 'table CSV first row');
    let bad = false;
    try { await V.decodeMap(Object.assign({}, m0, { nx: nx + 1 })); } catch (err) { bad = true; }
    check(bad, 'an image of the wrong size is an error');
  }

  // ---- repeat groups: the same statistics, tables and CSV as repeats.py ----
  if (fx.repeats) {
    const rx = fx.repeats;
    rx.stats.forEach((w, i) => {
      const got = V.repStats([[10.921, 10.251, 10.525], [5.0], [1.0, 3.0], []][i]);
      near(got.mean === null ? 0 : got.mean, w[0] === null ? 0 : w[0], 'stats mean ' + i);
      eq(got.sd === null ? null : +got.sd.toFixed(9), w[1] === null ? null : +w[1].toFixed(9), 'stats sd ' + i);
      eq(got.n, w[2], 'stats n ' + i);
    });
    rx.labels.forEach((c, i) => eq(V.repShortLabels(c[0]), c[1], 'short labels ' + i));
    ['plain', 'ticked'].forEach((which) => {
      const ex = rx[which];
      const include = {};
      ex.excluded.forEach((k) => { include[k] = false; });
      const [g] = V.repCompute(ex.pages, rx.specs, include, false, null, '');
      eq(g.members, ex.members, which + ': members');
      eq(g.notes, ex.notes, which + ': notes');
      eq(g.rows.length, ex.rows.length, which + ': row count');
      g.rows.forEach((r, i) => {
        const w = ex.rows[i];
        eq([r.region, r.state, r.n, r.of], [w.region, w.state, w.n, w.of], which + ': row ' + i + ' identity');
        near(r.mean, w.mean, which + ': mean ' + i);
        if (w.sd === null) eq(r.sd, null, which + ': no sd ' + i); else near(r.sd, w.sd, which + ': sd ' + i);
        eq(r.values.map((v) => v[0]), w.values.map((v) => v[0]), which + ': value labels ' + i);
        r.values.forEach((v, k) => near(v[1], w.values[k][1], which + ': value ' + i + ',' + k));
      });
      // the tables (a column per member, or none) and the CSV
      [['table6', 6], ['table2', 2]].forEach(([key, max]) => {
        const t = V.repTable(g, max), w = ex[key];
        eq(t.header, w[0], which + ': ' + key + ' header');
        eq(t.rows.map((r) => [r.kind, r.cells, r.partial]), w[1], which + ': ' + key + ' rows');
      });
      eq(V.repCsv([g]), ex.csv, which + ': CSV');
    });
    // the tick changes the numbers: O 1s of spot 2 is out, so the Ti 2p share rises there
    const a = V.repCompute(rx.plain.pages, rx.specs, {}, false, null, '')[0];
    const b = V.repCompute(rx.ticked.pages, rx.specs, { 's1r0:1': false }, false, null, '')[0];
    check(Math.abs(a.rows[0].mean - b.rows[0].mean) > 1e-6, 'unticking a region moves the group mean');
    // a group with one member says it has no spread, and the page shows only groups of two or more
    const one = V.repCompute(rx.plain.pages, [{ name: 'One', kind: 'locations', members: ['s0'] }], {}, false, null, '');
    check(one[0].notes.some((n) => n.indexOf('Only one member') === 0), 'one member says so');
    eq(V.repShown(one).length, 0, 'a one-member group is not shown');
    eq(V.repShown([a]).length, 1, 'a group of three is shown');
    // a sid the page does not have is a note, not an error
    const gone = V.repCompute(rx.plain.pages, [{ name: 'G', kind: 'nonsense', members: ['s0', 'zz'] }], {}, false, null, '')[0];
    eq(gone.kind, 'locations', 'an unknown kind reads as locations');
    check(gone.notes.some((n) => n.indexOf('not among') >= 0), 'a missing sample is noted');
    // Download everything carries the file, honouring the ticks
    const bundle = V.bundleFiles({ samples: [], details: {}, files: [], repeat_groups: rx.specs, generated: '2026-01-01T00:00:00' }, [], { include: {} });
    check(bundle.every((f) => f.name !== 'quantification_repeat_groups.csv'), 'no quantification, no repeat file');
  }

  // ---- the Columns / Shape boxes (the same shapes as the desktop's) ----
  check(V.gridColumns('auto') === '' && V.gridColumns('0') === '' && V.gridColumns('9') === '', 'auto columns leave the CSS alone');
  check(V.gridColumns('2') === 'repeat(2, minmax(0, 1fr))', 'two columns');
  check(V.PANEL_SHAPES.p34 === 4 / 3 && V.PANEL_SHAPES.p23 === 1.5 && V.PANEL_SHAPES.square === 1 && V.PANEL_SHAPES.l43 === 0.75 && V.PANEL_SHAPES.auto === null, 'panel shapes');

  // ---- hostile text stays text: nothing in the pure half builds HTML ----
  const src = fs.readFileSync(path.join(__dirname, '..', 'viewer', 'viewer.js'), 'utf8');
  check(src.indexOf('innerHTML') < 0, 'viewer never uses innerHTML');
  check(src.indexOf('insertAdjacentHTML') < 0, 'viewer never uses insertAdjacentHTML');
  check(src.indexOf('eval(') < 0 && src.indexOf('new Function') < 0, 'viewer never evals');
  check(!/https?:\/\//.test(src), 'viewer mentions no URLs');
  ['__DATA__', '__CSS__', '__JS__', '__TITLE__'].forEach((t) => check(src.indexOf(t) < 0, 'viewer.js does not contain ' + t));

  console.log((failed ? 'FAILED ' : 'ok ') + (ran - failed) + '/' + ran + ' checks');
  process.exit(failed ? 1 : 0);
})().catch((e) => { console.log('ERROR ' + (e && e.stack || e)); process.exit(2); });
