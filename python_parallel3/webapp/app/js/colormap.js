// colormap.js — palettes. viridis is the colorblind-safe default. Each entry maps
// t in [0,1] -> [r,g,b] (0-255) by interpolating control points.
window.VIZ = window.VIZ || {};
(function () {
  function ramp(stops) {
    return function (t) {
      t = Math.max(0, Math.min(1, t)); const s = t * (stops.length - 1);
      const i = Math.min(stops.length - 2, Math.floor(s)), f = s - i, a = stops[i], b = stops[i + 1];
      return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
    };
  }
  VIZ.colormaps = {
    viridis: ramp([[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]]),
    turbo:   ramp([[48, 18, 59], [70, 134, 251], [42, 232, 129], [249, 192, 41], [122, 4, 3]]),
    gray:    ramp([[30, 30, 34], [235, 235, 240]]),
    coolwarm: ramp([[59, 76, 192], [221, 221, 221], [180, 4, 38]])
  };
  VIZ.COLORBLIND_SAFE = ['viridis'];
})();
