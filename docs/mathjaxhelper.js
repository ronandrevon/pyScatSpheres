MathJax.Hub.Config({
  config: ["MMLorHTML.js"],
  jax: ["input/TeX", "output/HTML-CSS", "output/NativeMML"],
  extensions: ["MathMenu.js", "MathZoom.js"],
  TeX: {
    Macros: {
        cc: ["{\\mathcal #1}",1],
        braket: ["{\\langle #1 \\rangle}",1],
        bb: ["{\\mathbf #1}",1],
        dP: "{\\partial}",
        grad: "{\\nabla}",
        hn:"{h_n^{(1)}}",
        hl:"{h_l^{(1)}}",
        hbsm:"{\\frac{\\hbar^2}{2m_0}}",
        hbsme:["{\\frac{\\hbar^2}{2m_{#1}}}",1],
        vecTwo: ["{\\left[\\begin{array}{c} #1 \\\\ #2 \\end{array}\\right]}",2],
        matTwo: ["{\\left[\\begin{array}{cc} #1 & #2 \\\\ #3 & #4 \\end{array}\\right]}",4],
        matFour:["{\\left[\\begin{array}{cccc} #1\\\\ #2\\\\ #3\\\\ #4 \\end{array}\\right]}",4],
    }
  }
});
