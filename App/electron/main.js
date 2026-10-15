// Magibu - Electron shell: starts R (bundled) with the Shiny app on a local port
// and shows it in a window. Stops R on exit.
const { app, BrowserWindow, dialog, ipcMain, shell } = require('electron');
const { spawn, execFile } = require('child_process');
const path = require('path');
const fs = require('fs');
const net = require('net');
const http = require('http');

const packaged = app.isPackaged;
const root = packaged ? process.resourcesPath : path.join(__dirname, '..');
const appDir = path.join(root, 'app');
const modelsDir = path.join(root, 'models');
const fwd = (p) => p.replace(/\\/g, '/');

let rProc = null, win = null, quitting = false, logFile = null;

function findRscript() {
  if (process.env.MAGIBU_RSCRIPT) return process.env.MAGIBU_RSCRIPT;
  const bundled = path.join(root, packaged ? 'R' : path.join('build', 'R-portable'), 'bin',
                            process.platform === 'win32' ? 'Rscript.exe' : 'Rscript');
  return fs.existsSync(bundled) ? bundled : 'Rscript';   // no bundled R: use the one installed on the computer
}

function freePort() {
  return new Promise((resolve, reject) => {
    const s = net.createServer();
    s.once('error', reject);
    s.listen(0, '127.0.0.1', () => { const p = s.address().port; s.close(() => resolve(p)); });
  });
}

function waitForServer(url, timeoutMs) {
  const t0 = Date.now();
  return new Promise((resolve, reject) => {
    const tick = () => {
      if (!rProc) return reject(new Error('R stopped during start-up.'));
      http.get(url, (res) => { res.resume(); resolve(); }).on('error', () => {
        if (Date.now() - t0 > timeoutMs) reject(new Error('Timed out waiting for R to start.'));
        else setTimeout(tick, 400);
      });
    };
    tick();
  });
}

function logTail(n) {
  try { return fs.readFileSync(logFile, 'utf8').split(/\r?\n/).slice(-n).join('\n'); } catch (e) { return ''; }
}

function stopR() {
  if (!rProc) return;
  const pid = rProc.pid; rProc = null;
  if (process.platform === 'win32') execFile('taskkill', ['/pid', String(pid), '/T', '/F'], () => {});
  else { try { process.kill(-pid, 'SIGTERM'); } catch (e) { try { process.kill(pid, 'SIGTERM'); } catch (e2) {} } }
}

async function start() {
  logFile = path.join(app.getPath('userData'), 'magibu.log');
  fs.mkdirSync(path.dirname(logFile), { recursive: true });
  win = new BrowserWindow({
    width: 1500, height: 950, show: !process.env.MAGIBU_TEST_SHOT, autoHideMenuBar: true, title: 'Magibu',
    webPreferences: { preload: path.join(__dirname, 'preload.js'), contextIsolation: true, nodeIntegration: false }
  });
  win.on('closed', () => { win = null; });
  await win.loadFile(path.join(__dirname, 'splash.html'));

  const port = await freePort();
  const url = `http://127.0.0.1:${port}`;
  const out = fs.openSync(logFile, 'w');
  const env = Object.assign({}, process.env,
                            { MAGIBU_MODELS: modelsDir, MAGIBU_APP_DIR: fwd(appDir), MAGIBU_PORT: String(port) });
  if (packaged) { env.R_LIBS_USER = path.join(root, 'R', 'no-user-library'); delete env.R_LIBS; delete env.R_HOME; }
  rProc = spawn(findRscript(), ['--vanilla', path.join(appDir, 'launch.R')],
                { env, stdio: ['ignore', out, out], windowsHide: true, detached: process.platform !== 'win32' });
  rProc.on('error', (e) => { rProc = null; fail('Cannot start R: ' + e.message); });
  rProc.on('exit', (c) => { if (rProc && !quitting) { rProc = null; fail('R stopped (code ' + c + ').'); } });

  try {
    await waitForServer(url, 180000);
    if (!win) return;
    await win.loadURL(url);
    if (process.env.MAGIBU_TEST_SHOT) testShot();
  } catch (e) { fail(e.message); }
}

function fail(msg) {
  const detail = logTail(25);
  if (process.env.MAGIBU_TEST_SHOT) { console.error('FAIL', msg, '\n', detail); app.exit(1); return; }
  if (win) win.webContents.executeJavaScript(
    `document.getElementById('msg').textContent=${JSON.stringify(msg)};document.getElementById('err').textContent=${JSON.stringify(detail)};`
  ).catch(() => {});
  else dialog.showErrorBox('Magibu', msg + '\n\n' + detail);
}

// automatic test (no screen): checks the native bridge, saves a screenshot of the window, exits
async function testShot() {
  await new Promise((r) => setTimeout(r, 8000));
  const bridge = await win.webContents.executeJavaScript(
    "typeof window.magibuNative + '|' + document.body.classList.contains('mg-browser') + '|' + document.title");
  const img = await win.webContents.capturePage();
  fs.writeFileSync(process.env.MAGIBU_TEST_SHOT, img.toPNG());
  console.log('TEST bridge|browserMode|title =', bridge);
  quitting = true; stopR(); app.exit(0);
}

ipcMain.handle('magibu:pick', async (_e, mode) => {
  const props = mode === 'dir' ? ['openDirectory', 'createDirectory']
              : mode === 'files' ? ['openFile', 'multiSelections'] : ['openFile'];
  const r = await dialog.showOpenDialog(win, { properties: props });
  return r.canceled ? [] : r.filePaths.map(fwd);
});
ipcMain.handle('magibu:openPath', (_e, p) => shell.openPath(p));

if (!app.requestSingleInstanceLock()) app.quit();
else {
  app.on('second-instance', () => { if (win) { if (win.isMinimized()) win.restore(); win.focus(); } });
  app.whenReady().then(start);
  app.on('window-all-closed', () => { quitting = true; stopR(); app.quit(); });
  app.on('before-quit', () => { quitting = true; stopR(); });
}
