// Bridge between the interface and the native dialogs of the executable (window.magibuNative, exposed by Electron).
// In a normal browser the app works the same: paths are typed by hand.
(function () {
  var outDir = null;
  $(function () { if (!window.magibuNative) $('body').addClass('mg-browser'); });

  $(document).on('click', '.mg-browse', async function () {
    if (!window.magibuNative) return;
    var target = $(this).data('target'), mode = $(this).data('mode');
    if (mode === 'pred') mode = $('input[name=pred_mode]:checked').val();   // file | dir | files
    var picked = await window.magibuNative.pick(mode);
    if (!picked || picked.length === 0) return;
    $('#' + target).val(picked.join('\n')).trigger('change');
  });

  $(document).on('click', '.mg-open', function () {
    if (window.magibuNative && outDir) window.magibuNative.openPath(outDir);
  });

  $(document).on('shiny:connected', function () {
    Shiny.addCustomMessageHandler('mg-outdir', function (d) { outDir = d; });
  });
})();
