// Exposes only two native functions to the interface: choosing files/folders and opening a folder.
const { contextBridge, ipcRenderer } = require('electron');
contextBridge.exposeInMainWorld('magibuNative', {
  pick: (mode) => ipcRenderer.invoke('magibu:pick', mode),
  openPath: (p) => ipcRenderer.invoke('magibu:openPath', p)
});
