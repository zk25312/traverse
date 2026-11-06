/**
 * 导线控制点坐标平差 - 后端服务
 * 零依赖 Node.js：静态页面 + 数据共享存储 API
 * 数据存 NAS 共享目录，所有设备（手机/平板/电脑）读写同一份
 */
const http = require("http");
const fs = require("fs");
const path = require("path");
const PORT = process.env.PORT || 5001;
const SOCKET_PATH = process.env.SOCKET_PATH || "";
const GATEWAY_PREFIX = process.env.GATEWAY_PREFIX || "/app/traverse";
const DATA_DIR = process.env.DATA_DIR || path.join(__dirname, "data");
const DATA_FILE = path.join(DATA_DIR, "survey-data.json");
const STATIC_DIR = process.env.FRONTEND_DIST || path.join(__dirname, "public");
fs.mkdirSync(DATA_DIR, { recursive: true });
/* ---------- 数据读写（原子写入，防损坏） ---------- */
function loadData() {
  try {
    if (fs.existsSync(DATA_FILE)) {
      const d = JSON.parse(fs.readFileSync(DATA_FILE, "utf8"));
      if (d && typeof d === "object") return d;
    }
  } catch (e) { console.error("load failed:", e.message); }
  return { knowns: [], setups: [], version: 0, updatedAt: null };
}
function saveData(d) {
  d.version = (d.version || 0) + 1;
  d.updatedAt = new Date().toISOString();
  const tmp = DATA_FILE + ".tmp";
  fs.writeFileSync(tmp, JSON.stringify(d, null, 2), "utf8");
  fs.renameSync(tmp, DATA_FILE);
}
function sendJson(res, code, obj) {
  res.writeHead(code, {
    "Content-Type": "application/json; charset=utf-8",
    "Cache-Control": "no-store",
  });
  res.end(JSON.stringify(obj));
}
/* ---------- API ---------- */
function handleApi(req, res, pathname) {
  /* 读取共享数据 */
  if (req.method === "GET" && pathname === "/api/data") {
    return sendJson(res, 200, loadData());
  }
  /* 保存共享数据（整份覆盖） */
  if (req.method === "POST" && pathname === "/api/data") {
    let raw = "";
    req.on("data", (c) => {
      raw += c;
      if (raw.length > 8 * 1024 * 1024) { req.destroy(); }   // 上限 8MB
    });
    req.on("end", () => {
      let body;
      try { body = JSON.parse(raw || "{}"); } catch (e) { return sendJson(res, 400, { error: "invalid json" }); }
      if (!body || !Array.isArray(body.knowns) || !Array.isArray(body.setups)) {
        return sendJson(res, 400, { error: "invalid data structure" });
      }
      const d = { knowns: body.knowns, setups: body.setups, version: loadData().version || 0 };
      saveData(d);
      sendJson(res, 200, { ok: true, version: d.version, updatedAt: d.updatedAt });
    });
    return;
  }
  /* 清空共享数据 */
  if (req.method === "POST" && pathname === "/api/clear") {
    const d = { knowns: [], setups: [], version: 0 };
    saveData(d);
    return sendJson(res, 200, { ok: true });
  }
  sendJson(res, 404, { error: "not found" });
}
/* ---------- 静态文件 ---------- */
const MIME = {
  ".html": "text/html; charset=utf-8", ".css": "text/css; charset=utf-8",
  ".js": "application/javascript; charset=utf-8", ".json": "application/json; charset=utf-8",
  ".png": "image/png", ".jpg": "image/jpeg", ".svg": "image/svg+xml", ".ico": "image/x-icon",
};
function serveStatic(res, relPath) {
  let target = path.normalize(path.join(STATIC_DIR, relPath));
  if (!target.startsWith(path.normalize(STATIC_DIR))) {
    res.writeHead(400); return res.end("Bad Request");
  }
  if (relPath === "/" || relPath === "" || !fs.existsSync(target) || fs.statSync(target).isDirectory()) {
    if (path.extname(relPath) === "") {
      target = path.join(STATIC_DIR, "index.html");
    } else {
      res.writeHead(404, { "Content-Type": "text/plain; charset=utf-8" });
      return res.end("404 Not Found");
    }
  }
  if (!fs.existsSync(target)) {
    res.writeHead(404, { "Content-Type": "text/plain; charset=utf-8" });
    return res.end("404 Not Found");
  }
  const ext = path.extname(target).toLowerCase();
  res.writeHead(200, {
    "Content-Type": MIME[ext] || "application/octet-stream",
    "Cache-Control": "no-cache",
  });
  fs.createReadStream(target).pipe(res);
}
/* ---------- 路由 ---------- */
const server = http.createServer((req, res) => {
  let pathname = decodeURIComponent((req.url || "/").split("?")[0]);
  if (pathname.startsWith(GATEWAY_PREFIX + "/")) {
    pathname = pathname.slice(GATEWAY_PREFIX.length);
  } else if (pathname === GATEWAY_PREFIX) {
    res.writeHead(301, { Location: GATEWAY_PREFIX + "/" });
    return res.end();
  }
  if (pathname.startsWith("/api/")) return handleApi(req, res, pathname);
  return serveStatic(res, pathname);
});
if (SOCKET_PATH) {
  try { fs.rmSync(SOCKET_PATH, { force: true }); } catch (e) { /* ignore */ }
  server.listen(SOCKET_PATH, () => console.log(`traverse running on socket ${SOCKET_PATH}`));
} else {
  server.listen(PORT, () => console.log(`traverse running at http://localhost:${PORT}`));
}