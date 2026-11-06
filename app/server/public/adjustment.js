/*
 * 导线控制点坐标平差核心算法（纯计算，浏览器/Node 通用）
 * 外业记录：每站观测各目标的 水平角读数(即方位角读数) + 水平距离
 * 平差：方向(方位角)+边长 间接平差（未知数：待定点X/Y + 各测站定向角）
 * 高程：由已知点高差按距离比例推算（无竖直角时）
 * 支持：闭合导线 / 附合导线 / 支导线 / 往返支导线
 */
(function (root) {
  "use strict";
  var D2R = Math.PI / 180;
  var R2D = 180 / Math.PI;
  var RO = 206264.80624709636;
  function isNum(v) { return v !== "" && v !== null && v !== undefined && isFinite(+v); }
  function normAng(a) { while (a < 0) a += 2 * Math.PI; while (a >= 2 * Math.PI) a -= 2 * Math.PI; return a; }
  function angDiff(a, b) { var d = a - b; while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI; return d; }
  /* 角度解析：ddd.mmss（测量惯例）/ "65°47′34.1″" / "65 47 34.1" / 十进制度 */
  function parseDMS(str) {
    str = String(str == null ? "" : str).trim();
    if (!str) return NaN;
    var m = str.match(/^(-?\d+(?:\.\d+)?)\s*[°度]\s*(\d+(?:\.\d+)?)\s*['′分]\s*(\d+(?:\.\d+)?)\s*["″秒]?$/);
    if (m) return (+m[1]) + (+m[2]) / 60 + (+m[3]) / 3600;
    m = str.match(/^(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)$/);
    if (m) return (+m[1]) + (+m[2]) / 60 + (+m[3]) / 3600;
    var v = parseFloat(str);
    if (isNaN(v)) return NaN;
    var neg = v < 0; v = Math.abs(v);
    var d = Math.floor(v);
    var fr = (v - d) * 100;
    var mi = Math.floor(fr + 1e-9);
    var ss = (fr - mi) * 100;
    if (mi >= 60 || ss >= 60.001) return NaN;
    var deg = d + mi / 60 + ss / 3600;
    return neg ? -deg : deg;
  }
  function fmtDMS(deg, sd) {
    if (!isFinite(deg)) return "";
    sd = sd == null ? 1 : sd;
    var neg = deg < 0; deg = Math.abs(deg);
    var d = Math.floor(deg);
    var mf = (deg - d) * 60;
    var mi = Math.floor(mf + 1e-9);
    var ss = (mf - mi) * 60;
    var digits = Math.round((mi * 100 + ss) * Math.pow(10, sd));
    var frac = String(digits);
    while (frac.length < 4 + sd) frac = "0" + frac;
    return (neg ? "-" : "") + d + "." + frac;
  }
  function solveLinear(A, b) {
    var n = b.length;
    var M = A.map(function (r) { return r.slice(); });
    var x = b.slice();
    for (var c = 0; c < n; c++) {
      var p = c;
      for (var r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      if (Math.abs(M[p][c]) < 1e-14) return null;
      if (p !== c) { var tm = M[p]; M[p] = M[c]; M[c] = tm; var tb = x[p]; x[p] = x[c]; x[c] = tb; }
      var pv = M[c][c];
      for (var r2 = c + 1; r2 < n; r2++) {
        var f = M[r2][c] / pv;
        if (!f) continue;
        for (var cc = c; cc < n; cc++) M[r2][cc] -= f * M[c][cc];
        x[r2] -= f * x[c];
      }
    }
    for (var i = n - 1; i >= 0; i--) {
      var s = x[i];
      for (var j = i + 1; j < n; j++) s -= M[i][j] * x[j];
      x[i] = s / M[i][i];
    }
    return x;
  }
  function invertMatrix(A) {
    var n = A.length, cols = [];
    for (var j = 0; j < n; j++) {
      var e = new Array(n); e.fill(0); e[j] = 1;
      var x = solveLinear(A, e);
      if (!x) return null;
      cols.push(x);
    }
    var Q = [];
    for (var i = 0; i < n; i++) Q.push(new Array(n));
    for (i = 0; i < n; i++) for (var j2 = 0; j2 < n; j2++) Q[i][j2] = cols[j2][i];
    return Q;
  }
  function solveLS(rows, t) {
    if (t <= 0 || !rows.length) return null;
    var N = [], rhs = new Array(t); rhs.fill(0);
    for (var i = 0; i < t; i++) { N.push(new Array(t)); N[i].fill(0); }
    rows.forEach(function (r) {
      var ks = Object.keys(r.coef);
      for (var a = 0; a < ks.length; a++) {
        var ia = +ks[a], ca = r.coef[ks[a]];
        rhs[ia] += r.w * ca * r.l;
        for (var b = 0; b < ks.length; b++) N[ia][+ks[b]] += r.w * ca * r.coef[ks[b]];
      }
    });
    var x = solveLinear(N, rhs);
    if (!x) return null;
    var Q = invertMatrix(N);
    if (!Q) return null;
    var vtPv = 0;
    rows.forEach(function (r) {
      var s = 0;
      for (var k in r.coef) s += r.coef[k] * x[+k];
      r.v = s - r.l;
      vtPv += r.w * r.v * r.v;
    });
    return { x: x, Q: Q, vtPv: vtPv };
  }
  function resectionTwoDist(p1, r1, p2, r2, rays) {
    var dx = p2.x - p1.x, dy = p2.y - p1.y;
    var dd = Math.hypot(dx, dy);
    if (dd < 1e-9) return null;
    if (r1 + r2 < dd - 1e-9 || Math.abs(r1 - r2) > dd + 1e-9) return null;
    var a = (dd * dd + r1 * r1 - r2 * r2) / (2 * dd);
    var h2 = r1 * r1 - a * a;
    var h = h2 > 0 ? Math.sqrt(h2) : 0;
    var ux = dx / dd, uy = dy / dd;
    var bx = p1.x + a * ux, by = p1.y + a * uy;
    var cands = [{ x: bx - h * uy, y: by + h * ux }];
    if (h > 1e-9) cands.push({ x: bx + h * uy, y: by - h * ux });
    if (cands.length === 1) return cands[0];
    function fit(c) {
      var s = 0, n = 0;
      for (var i = 0; i < rays.length; i++) for (var j = i + 1; j < rays.length; j++) {
        var A1 = Math.atan2(rays[i].py - c.y, rays[i].px - c.x);
        var A2 = Math.atan2(rays[j].py - c.y, rays[j].px - c.x);
        s += Math.pow(angDiff(A2 - A1, rays[j].dir - rays[i].dir), 2); n++;
      }
      return n ? s / n : 0;
    }
    return fit(cands[0]) <= fit(cands[1]) ? cands[0] : cands[1];
  }
  function lineIntersect(p1, b1, p2, b2) {
    var u1x = Math.cos(b1), u1y = Math.sin(b1);
    var u2x = Math.cos(b2), u2y = Math.sin(b2);
    var cross = u1x * u2y - u1y * u2x;
    if (Math.abs(cross) < 1e-12) return null;
    var dx = p2.x - p1.x, dy = p2.y - p1.y;
    var t = (dx * u2y - dy * u2x) / cross;
    return { x: p1.x + u1x * t, y: p1.y + u1y * t };
  }
  /* ================= 主入口 ================= */
  function computeAdjustment(input) {
    var errors = [], warnings = [];
    var cfg = input.config || {};
    var saSec = isNum(cfg.saSec) ? +cfg.saSec : 2.5;
    var aMM = isNum(cfg.aMM) ? +cfg.aMM : 2;
    var bPPM = isNum(cfg.bPPM) ? +cfg.bPPM : 2;
    var dirDecimal = !!cfg.dirDecimal;
    var sDirRad = saSec / RO;
    function parseDir(v) {
      if (v === "" || v == null || isNaN(+v)) return NaN;
      if (typeof v === "number") return v;
      return dirDecimal ? parseFloat(v) : parseDMS(v);
    }
    var pts = {};
    function getPt(n) { n = String(n).trim(); if (!pts[n]) pts[n] = { name: n }; return pts[n]; }
    (input.knowns || []).forEach(function (k) {
      if (!k || !String(k.name || "").trim()) return;
      var p = getPt(k.name);
      if (isNum(k.X) && isNum(k.Y)) { p.X = +k.X; p.Y = +k.Y; p.knownXY = true; }
      else if (isNum(k.X) || isNum(k.Y)) errors.push("已知点 " + p.name + "：X 与 Y 必须同时填写");
      if (isNum(k.H)) { p.H = +k.H; p.knownH = true; }
    });
    /* ---- 测站与观测 ---- */
    var setups = [], setupMap = {};
    (input.setups || []).forEach(function (s) {
      var sn = String(s.station || "").trim();
      if (!sn) return;
      var st = setupMap[sn];
      if (!st) { st = setupMap[sn] = { station: sn, obs: [], relH: null }; setups.push(st); }
      if (st.relH == null && isNum(s.relH)) st.relH = +s.relH;
      (s.obs || []).forEach(function (o) {
        var tn = String(o.target || "").trim();
        if (!tn || tn === sn) return;
        var dir = isNum(o.dir) ? normAng(parseDir(o.dir) * D2R) : null;
        var D = isNum(o.dist) ? +o.dist : null;
        if (D != null && D <= 0) { warnings.push("测站 " + sn + " → " + tn + "：距离非正，已忽略"); D = null; }
        var relH = isNum(o.relH) ? +o.relH : null;
        st.obs.push({ target: getPt(tn), dir: dir, D: D, face: o.face || "L", relH: relH });
      });
    });
    setups = setups.filter(function (st) { return st.obs.length > 0; });
    /* ---- 正倒镜一测回合并：同一目标的正镜/倒镜取平均（倒镜 ±180°） ---- */
    setups.forEach(function (st) {
      var byTarget = {}, order = [];
      st.obs.forEach(function (o) {
        var key = o.target.name;
        if (!byTarget[key]) { byTarget[key] = { target: o.target, L: null, R: null, D: null, Ds: [], relHs: [] }; order.push(key); }
        var g = byTarget[key];
        if (o.dir != null) {
          if (o.face === "R") g.R = o.dir; else g.L = o.dir;
        }
        if (o.D != null) g.Ds.push(o.D);
        if (o.relH != null) g.relHs.push(o.relH);
      });
      var merged = [];
      order.forEach(function (key) {
        var g = byTarget[key];
        var dir = null;
        if (g.L != null && g.R != null) {
          /* 一测回平均：倒镜读数归算到正镜方向（差接近±180°时±180°） */
          var rAdj = g.R;
          var d = angDiff(rAdj, g.L);
          if (Math.abs(Math.abs(d) - Math.PI) < Math.PI / 2) {
            /* 让倒镜靠近正镜：差为正则减180，为负则加180 */
            rAdj = d > 0 ? rAdj - Math.PI : rAdj + Math.PI;
          }
          /* 关键：先归一化到 [0,2π) 再平均，避免跨 0/360 边界 */
          rAdj = normAng(rAdj);
          var Ln = normAng(g.L);
          /* 若两者仍相差接近 2π（跨零），把小的加上 2π 再平均 */
          if (Math.abs(Ln - rAdj) > Math.PI) {
            if (Ln < rAdj) Ln += 2 * Math.PI; else rAdj += 2 * Math.PI;
          }
          dir = normAng((Ln + rAdj) / 2);
        } else if (g.L != null) dir = g.L;
        else if (g.R != null) dir = g.R;
        var D = g.Ds.length ? g.Ds.reduce(function (a, b) { return a + b; }, 0) / g.Ds.length : null;
        var relH = g.relHs.length ? g.relHs[0] : null;
        merged.push({ target: g.target, dir: dir, D: D, relH: relH });
      });
      st.obs = merged;
    });
    if (!setups.length) { errors.push("没有有效的观测记录"); return fail(); }
    var nKnownXY = Object.keys(pts).filter(function (n) { return pts[n].knownXY; }).length;
    if (nKnownXY < 2) { errors.push("已知平面点（XY）少于 2 个，无法确定坐标基准"); return fail(); }
    function hasXY(p) { return p.knownXY || isNum(p.X); }
    function hasH(p) { return p.knownH || isNum(p.H); }
    /* ---- 导线类型识别 ---- */
    var stationNames = setups.map(function (s) { return s.station; });
    var knownStations = stationNames.filter(function (n) { return getPt(n).knownXY; });
    var unknownStations = stationNames.filter(function (n) { return !getPt(n).knownXY; });
    var allTargets = {};
    setups.forEach(function (st) { st.obs.forEach(function (o) { allTargets[o.target.name] = o.target; }); });
    var knownTargets = Object.keys(allTargets).filter(function (n) { return allTargets[n].knownXY; });
    var traverseType = "未知";
    /* 统计：已知测站数、导线末端是否连到已知点 */
    var hasKnownTarget = setups.some(function (s) {
      return s.obs.some(function (o) { return o.target.knownXY; });
    });
    if (knownStations.length >= 2) {
      var first = knownStations[0], last = knownStations[knownStations.length - 1];
      var isClosed = (first === last) || setups.some(function (s) {
        return s.station === first && s.obs.some(function (o) { return o.target.name === last; });
      });
      traverseType = isClosed ? "闭合导线" : "附合导线";
    } else if (knownStations.length === 1) {
      /* 1个已知测站：若导线末端连到已知点 → 附合导线；若未知站回测已知站 → 往返支导线；否则支导线 */
      if (hasKnownTarget) {
        traverseType = "附合导线";
      } else {
        var back = setups.some(function (s) {
          return unknownStations.indexOf(s.station) >= 0 &&
            s.obs.some(function (o) { return o.target.knownXY; });
        });
        traverseType = back ? "往返支导线" : "支导线";
      }
    } else {
      traverseType = hasKnownTarget ? "附合导线" : "支导线";
    }
    /* ---- 近似坐标推算 ---- */
    /* 迭代：先用已知点定位能定位的站，再用已定位站传递/交会，直到不再变化 */
    function computeOrientation(st) {
      var sp = getPt(st.station);
      if (!hasXY(sp)) return false;
      var sum = 0, n = 0;
      for (var k = 0; k < st.obs.length; k++) {
        var o = st.obs[k];
        if (o.dir == null) continue;
        if (!hasXY(o.target)) continue;
        sum += angDiff(Math.atan2(o.target.Y - sp.Y, o.target.X - sp.X), o.dir);
        n++;
      }
      if (n > 0) { st.O0 = sum / n; return true; }
      return false;
    }
    function propagateTraverse() {
      var changed = false;
      /* 第一步：所有已定位站算定向角 */
      for (var a = 0; a < setups.length; a++) {
        var st = setups[a];
        if (st.O0 === undefined || st.O0 === null) {
          if (computeOrientation(st)) changed = true;
        }
      }
      /* 第二步：未定位站由已定位站（有定向角）传递坐标 */
      for (var b = 0; b < setups.length; b++) {
        var st2 = setups[b];
        var sp2 = getPt(st2.station);
        if (hasXY(sp2)) continue;
        for (var c = 0; c < setups.length; c++) {
          if (c === b) continue;
          var other = setups[c];
          var op = getPt(other.station);
          if (!hasXY(op)) continue;
          if (other.O0 === undefined || other.O0 === null) continue;
          for (var d = 0; d < other.obs.length; d++) {
            var o2 = other.obs[d];
            if (o2.target.name !== st2.station) continue;
            if (o2.D == null || o2.dir == null) continue;
            var A = other.O0 + o2.dir;
            sp2.X = op.X + o2.D * Math.cos(A);
            sp2.Y = op.Y + o2.D * Math.sin(A);
            changed = true;
            break;
          }
          if (hasXY(sp2)) break;
        }
      }
      return changed;
    }
    for (var iter = 0; iter < 100; iter++) {
      var changed = false;
      if (propagateTraverse()) changed = true;
      setups.forEach(function (st) {
        var sp = getPt(st.station);
        /* 未知站：用 2 个已定位目标做距离交会 */
        if (!hasXY(sp)) {
          var kn = st.obs.filter(function (o) { return o.D != null && hasXY(o.target); });
          if (kn.length >= 2) {
            var b1 = kn[0], b2 = kn[1], bestD = -1;
            for (var i = 0; i < kn.length; i++) for (var j = i + 1; j < kn.length; j++) {
              var dd = Math.hypot(kn[i].target.X - kn[j].target.X, kn[i].target.Y - kn[j].target.Y);
              if (dd > bestD) { bestD = dd; b1 = kn[i]; b2 = kn[j]; }
            }
            var rays = kn.map(function (o) { return { dir: o.dir, px: o.target.X, py: o.target.Y }; });
            var cc = resectionTwoDist({ x: b1.target.X, y: b1.target.Y }, b1.D,
                                      { x: b2.target.X, y: b2.target.Y }, b2.D, rays);
            if (cc) { sp.X = cc.x; sp.Y = cc.y; changed = true; }
          }
        }
        /* 已定位站：算定向角 + 极坐标推算未定位目标 */
        if (hasXY(sp)) {
          if (computeOrientation(st)) changed = true;
          if (st.O0 !== undefined && st.O0 !== null) {
            st.obs.forEach(function (o) {
              if (o.D != null && o.dir != null && !hasXY(o.target)) {
                var A = st.O0 + o.dir;
                o.target.X = sp.X + o.D * Math.cos(A);
                o.target.Y = sp.Y + o.D * Math.sin(A);
                changed = true;
              }
            });
          }
        }
      });
      setups.forEach(function (st) {
        var sp = getPt(st.station);
        if (!hasXY(sp) || st.O0 == null) return;
        st.obs.forEach(function (o) {
          if (o.dir == null || hasXY(o.target)) return;
          var b1 = st.O0 + o.dir;
          setups.forEach(function (st2) {
            var sp2 = getPt(st2.station);
            if (st2 === st || !hasXY(sp2) || st2.O0 == null) return;
            for (var q = 0; q < st2.obs.length; q++) {
              var o2 = st2.obs[q];
              if (o2.target !== o.target || o2.dir == null || hasXY(o2.target)) continue;
              var ip = lineIntersect({ x: sp.X, y: sp.Y }, b1, { x: sp2.X, y: sp2.Y }, st2.O0 + o2.dir);
              if (ip) { o.target.X = ip.x; o.target.Y = ip.y; changed = true; }
              break;
            }
          });
        });
      });
      if (!changed) break;
    }
    var referenced = {};
    setups.forEach(function (st) {
      referenced[st.station] = getPt(st.station);
      st.obs.forEach(function (o) { referenced[o.target.name] = o.target; });
    });
    var missingXY = Object.keys(referenced).filter(function (n) { return !hasXY(referenced[n]); });
    if (missingXY.length) {
      errors.push("以下点无法求得近似坐标：" + missingXY.join("、") +
        "。要求：未知测站至少观测 2 个已知(或已算出)点且含距离；每个待定点至少有 1 条含距离的观测，或被 2 个测站的方向交会。");
      return fail();
    }
    /* ---- 平面间接平差 ---- */
    var uPts = [], uIdx = {};
    function addU(n) { if (uIdx[n] == null) { uIdx[n] = uPts.length; uPts.push(n); } }
    setups.forEach(function (st) {
      if (!getPt(st.station).knownXY) addU(st.station);
      st.obs.forEach(function (o) { if (!o.target.knownXY) addU(o.target.name); });
    });
    var nU2 = uPts.length * 2;
    var oriSetups = setups.filter(function (st) { return st.obs.some(function (o) { return o.dir != null; }); });
    var tPlan = nU2 + oriSetups.length;
    var dirRows = [], distRows = [];
    setups.forEach(function (st) {
      var sp = getPt(st.station);
      var sum = 0, n = 0;
      st.obs.forEach(function (o) {
        if (o.dir != null && hasXY(o.target)) { sum += angDiff(Math.atan2(o.target.Y - sp.Y, o.target.X - sp.X), o.dir); n++; }
      });
      st.O0 = n ? sum / n : 0;
      var oIdx = oriSetups.indexOf(st);
      st.obs.forEach(function (o) {
        var tp = o.target;
        var dx = tp.X - sp.X, dy = tp.Y - sp.Y;
        var D0 = Math.hypot(dx, dy);
        if (D0 <= 0) return;
        var A0 = Math.atan2(dy, dx);
        if (o.dir != null) {
          var coef = {};
          if (!tp.knownXY) { coef[2 * uIdx[tp.name]] = -Math.sin(A0) / D0; coef[2 * uIdx[tp.name] + 1] = Math.cos(A0) / D0; }
          if (!sp.knownXY) { coef[2 * uIdx[sp.name]] = Math.sin(A0) / D0; coef[2 * uIdx[sp.name] + 1] = -Math.cos(A0) / D0; }
          if (oIdx >= 0) coef[nU2 + oIdx] = -1;
          dirRows.push({ coef: coef, l: angDiff(o.dir, A0 - st.O0), w: 1,
            from: st.station, to: tp.name, obsDeg: o.dir * R2D, kind: "方向" });
        }
        if (o.D != null) {
          var coef2 = {};
          if (!tp.knownXY) { coef2[2 * uIdx[tp.name]] = Math.cos(A0); coef2[2 * uIdx[tp.name] + 1] = Math.sin(A0); }
          if (!sp.knownXY) { coef2[2 * uIdx[sp.name]] = -Math.cos(A0); coef2[2 * uIdx[sp.name] + 1] = -Math.sin(A0); }
          if (Object.keys(coef2).length) {
            var sig = (aMM + bPPM * D0 / 1000) / 1000;
            distRows.push({ coef: coef2, l: o.D - D0, w: (sDirRad / sig) * (sDirRad / sig),
              from: st.station, to: tp.name, obsM: o.D, kind: "边长" });
          }
        }
      });
    });
    var rowsPlan = dirRows.concat(distRows);
    var nPlan = rowsPlan.length;
    var resPlan = solveLS(rowsPlan, tPlan);
    if (!resPlan) {
      errors.push("平面法方程奇异：起算/观测数据不足或图形结构不佳（未知数 " + tPlan + " 个，观测 " + nPlan + " 个）");
      return fail();
    }
    var rPlan = nPlan - tPlan, s0sec;
    if (rPlan > 0) s0sec = Math.sqrt(resPlan.vtPv / rPlan) * RO;
    else { s0sec = saSec; warnings.push("平面多余观测 r = " + rPlan + " ≤ 0，单位权中误差按先验值 " + saSec + "″ 给出，仅供参考"); }
    uPts.forEach(function (n, i) { var p = getPt(n); p.X += resPlan.x[2 * i]; p.Y += resPlan.x[2 * i + 1]; });
    oriSetups.forEach(function (st, k) { st.O0 += resPlan.x[nU2 + k]; });
    var mErr = {};
    uPts.forEach(function (n, i) {
      mErr[n] = {
        mX: (s0sec / RO) * Math.sqrt(resPlan.Q[2 * i][2 * i]),
        mY: (s0sec / RO) * Math.sqrt(resPlan.Q[2 * i + 1][2 * i + 1])
      };
    });
    /* ---- 高程 ---- */
    var knownHs = Object.keys(referenced).filter(function (n) { return referenced[n].knownH; });
    /* 方法1（优先）：相对高程锚定 —— 已知高程点也填了相对高程时，
       基准差 = 已知H − 相对高程；其余点 绝对H = 相对高程 + 基准差 */
    var relHMap = {};
    setups.forEach(function (st) {
      if (st.relH != null && relHMap[st.station] == null) relHMap[st.station] = st.relH;
      st.obs.forEach(function (o) {
        if (o.relH != null && relHMap[o.target.name] == null) relHMap[o.target.name] = o.relH;
      });
    });
    var relOffset = null, relOffsetN = 0;
    knownHs.forEach(function (n) {
      if (relHMap[n] != null) {
        relOffset = (relOffset == null ? 0 : relOffset) + (referenced[n].H - relHMap[n]);
        relOffsetN++;
      }
    });
    if (relOffsetN > 0) {
      relOffset = relOffset / relOffsetN;
      var relCount = 0;
      Object.keys(relHMap).forEach(function (n) {
        var p = referenced[n];
        if (!p || p.knownH) return;
        p.H = relHMap[n] + relOffset;
        p.hFromRel = true;
        relCount++;
      });
      warnings.push("已按相对高程反算 " + relCount + " 个点的高程（基准差 " + relOffset.toFixed(4) + " m，锚定自 " + relOffsetN + " 个已知高程点）");
    }
    /* 方法2（兜底）：按已知点高差与距离比例推算 */
    if (knownHs.length >= 1) {
      var knownList = knownHs.map(function (n) { return referenced[n]; });
      /* 建立观测图（含距离），用最短路径把高程从已知点传到未知点 */
      var adj = {};
      function addEdge(a, b, dist) {
        (adj[a] = adj[a] || []).push({ to: b, d: dist });
        (adj[b] = adj[b] || []).push({ to: a, d: dist });
      }
      setups.forEach(function (st) {
        st.obs.forEach(function (o) { if (o.D != null) addEdge(st.station, o.target.name, o.D); });
      });
      /* 用已知点间的高差/距离得到平均坡度，再按距离推算 */
      var slope = 0, slopeN = 0;
      for (var i1 = 0; i1 < knownList.length; i1++) for (var j1 = i1 + 1; j1 < knownList.length; j1++) {
        var A = knownList[i1], B = knownList[j1];
        var dAB = Math.hypot(B.X - A.X, B.Y - A.Y);
        if (dAB > 1e-6) { slope += (B.H - A.H) / dAB; slopeN++; }
      }
      slope = slopeN ? slope / slopeN : 0;
      /* 逐点推算：以所有已有高程的点（含相对高程反算的）为源，BFS 按距离比例 */
      var sources = Object.keys(referenced).filter(function (n) { return hasH(referenced[n]); });
      var queue = sources.slice(), visited = {};
      sources.forEach(function (n) { visited[n] = true; });
      while (queue.length) {
        var cur = queue.shift();
        var cp = referenced[cur];
        (adj[cur] || []).forEach(function (e) {
          if (visited[e.to]) return;
          var tp = referenced[e.to];
          if (!tp) return;
          if (!tp.knownH && !hasH(tp)) {
            tp.H = cp.H + slope * e.d;
            tp.hApprox = true;
          }
          visited[e.to] = true;
          queue.push(e.to);
        });
      }
      var noH = Object.keys(referenced).filter(function (n) { return !hasH(referenced[n]); });
      if (noH.length) warnings.push("以下点高程无法推算（与已知高程点不连通）：" + noH.join("、"));
      if (slopeN === 0) warnings.push("仅 1 个已知高程点，无法确定坡度，未知点高程按该点高程赋值（仅供参考）");
    } else {
      warnings.push("未提供已知高程点，高程不参与计算");
    }
    /* ---- 闭合差统计：导线边坐标增量之和 vs 理论值 ---- */
    var closure = null;
    if (traverseType === "闭合导线" || traverseType === "附合导线") {
      var sumDx = 0, sumDy = 0, sumD = 0, nEdge = 0;
      setups.forEach(function (st) {
        st.obs.forEach(function (o) {
          if (o.D == null || o.dir == null) return;
          if (!setupMap[o.target.name]) return;   /* 只统计站→站的导线边 */
          var A = st.O0 + o.dir;
          sumDx += o.D * Math.cos(A);
          sumDy += o.D * Math.sin(A);
          sumD += o.D;
          nEdge++;
        });
      });
      /* 理论值：闭合导线为 0；附合导线为 末已知点 - 首已知点 */
      var theoDx = 0, theoDy = 0;
      if (traverseType === "附合导线") {
        /* 找导线首站（已知）与末端已知点 */
        var startName = null, endName = null;
        for (var si = 0; si < setups.length; si++) {
          if (getPt(setups[si].station).knownXY) { startName = setups[si].station; break; }
        }
        for (var ei = setups.length - 1; ei >= 0; ei--) {
          var eobs = setups[ei].obs;
          for (var ej = 0; ej < eobs.length; ej++) {
            if (eobs[ej].target.knownXY && setupMap[eobs[ej].target.name]) { endName = eobs[ej].target.name; break; }
          }
          if (endName) break;
        }
        if (!endName && knownStations.length >= 2) endName = knownStations[knownStations.length - 1];
        if (startName && endName) {
          var p0 = getPt(startName), p1 = getPt(endName);
          if (p0.knownXY && p1.knownXY) { theoDx = p1.X - p0.X; theoDy = p1.Y - p0.Y; }
        }
      }
      var fx = sumDx - theoDx, fy = sumDy - theoDy;
      closure = { sumDx: sumDx, sumDy: sumDy, sumD: sumD, nEdge: nEdge,
        theoDx: theoDx, theoDy: theoDy, fx: fx, fy: fy,
        fL: Math.hypot(fx, fy),
        rel: sumD > 0 ? Math.hypot(fx, fy) / sumD : null };
    }
    /* ---- 输出 ---- */
    var outPts = [];
    Object.keys(referenced).forEach(function (n) {
      var p = referenced[n];
      var m = mErr[n];
      outPts.push({
        name: n,
        X: hasXY(p) ? p.X : null, Y: hasXY(p) ? p.Y : null, H: hasH(p) ? p.H : null,
        knownXY: !!p.knownXY, knownH: !!p.knownH, hApprox: !!p.hApprox, hFromRel: !!p.hFromRel,
        mX: m ? m.mX * 1000 : null, mY: m ? m.mY * 1000 : null,
        mP: m ? Math.hypot(m.mX, m.mY) * 1000 : null,
        mH: null
      });
    });
    return {
      ok: true, errors: errors, warnings: warnings, traverseType: traverseType,
      points: outPts,
      plan: {
        n: nPlan, t: tPlan, r: rPlan, s0: s0sec, nDir: dirRows.length, nDist: distRows.length,
        dirRows: dirRows.map(function (r) { return { from: r.from, to: r.to, obs: r.obsDeg, v: r.v * RO }; }),
        distRows: distRows.map(function (r) { return { from: r.from, to: r.to, obs: r.obsM, v: r.v * 1000 }; })
      },
      height: { method: "按已知点高差与距离比例推算", knownCount: (knownHs || []).length },
      closure: closure,
      orientations: oriSetups.map(function (st) { return { station: st.station, O0deg: normAng(st.O0) * R2D }; })
    };
    function fail() {
      var out = [];
      Object.keys(referenced || {}).forEach(function (n) {
        var p = referenced[n];
        out.push({ name: n, X: hasXY(p) ? p.X : null, Y: hasXY(p) ? p.Y : null, H: hasH(p) ? p.H : null,
          knownXY: !!p.knownXY, knownH: !!p.knownH, mX: null, mY: null, mP: null, mH: null });
      });
      return { ok: false, errors: errors, warnings: warnings, traverseType: traverseType, points: out,
        plan: null, height: null, closure: null, orientations: [] };
    }
  }
  var api = { computeAdjustment: computeAdjustment, parseDMS: parseDMS, fmtDMS: fmtDMS,
    D2R: D2R, R2D: R2D, RO: RO };
  root.SurveyAdjust = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);