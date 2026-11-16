// Export PNG (scale 3) and SVG from the standalone Plotly HTML figures using headless Chromium.
// Visualization only: it loads the HTML files written by make_tubitak_3d.py and calls Plotly.toImage.
// Note: for 3D (WebGL) scenes Plotly's SVG export embeds the rendered scene as a raster image;
// only titles, colorbars and annotations are vector elements.
// Usage (repository root): node visualizations/tubitak/export_static.js
const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');

const DIR = __dirname;
const FIGS = [
  ['fig1_duzensizlik_potansiyeli_3b', 900, 760],
  ['fig2_pinn_dst_karsilastirma_3b', 1500, 760],
  ['fig3_pinn_dst_mutlak_fark_3b', 900, 760],
];

(async () => {
  const browser = await chromium.launch({
    args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'],
  });
  for (const [name, w, h] of FIGS) {
    const page = await browser.newPage({ viewport: { width: w + 40, height: h + 40 } });
    await page.goto('file://' + path.join(DIR, name + '.html'));
    await page.waitForFunction(() => document.querySelector('.plotly-graph-div') &&
                                     document.querySelector('.plotly-graph-div')._fullLayout);
    await page.waitForTimeout(1500);
    for (const fmt of ['png', 'svg']) {
      const url = await page.evaluate(async ([fmt, w, h]) => {
        const gd = document.querySelector('.plotly-graph-div');
        return await Plotly.toImage(gd, { format: fmt, width: w, height: h, scale: fmt === 'png' ? 3 : 1 });
      }, [fmt, w, h]);
      let buf;
      if (url.startsWith('data:image/svg+xml,')) buf = Buffer.from(decodeURIComponent(url.split(',')[1]));
      else buf = Buffer.from(url.split(',')[1], 'base64');
      fs.writeFileSync(path.join(DIR, `${name}.${fmt}`), buf);
      console.log('wrote', `${name}.${fmt}`, buf.length);
    }
    await page.close();
  }
  // preview of the improved animation (first frame), for documentation only
  const page = await browser.newPage({ viewport: { width: 1540, height: 860 } });
  await page.goto('file://' + path.join(DIR, 'wavefunction_animation_white.html'));
  await page.waitForFunction(() => document.querySelector('.plotly-graph-div')?._fullLayout);
  await page.waitForTimeout(1500);
  await page.evaluate(() => Plotly.animate(document.querySelector('.plotly-graph-div'), ['k30'],
                                           { mode: 'immediate', frame: { duration: 0, redraw: true } }));
  await page.waitForTimeout(1500);
  await page.screenshot({ path: path.join(DIR, 'preview_wavefunction_animation_white_t_end.png') });
  await page.evaluate(() => { const gd = document.querySelector('.plotly-graph-div'); const b = gd._fullLayout.updatemenus[1].buttons[1]; return Plotly.update(gd, b.args[0], b.args[1]); });
  await page.waitForTimeout(1500);
  await page.screenshot({ path: path.join(DIR, 'preview_wavefunction_animation_white_t_end_closeup.png') });
  console.log('wrote animation preview');
  await browser.close();
})();
