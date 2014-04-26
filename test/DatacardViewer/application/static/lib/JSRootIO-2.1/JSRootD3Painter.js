// JSROOTD3Painter.js
//
// core methods for Javascript ROOT Graphics, using d3.js.
//

// The "source_dir" variable is defined in JSRootInterface.js

var d_tree, key_tree, defs;

var kWhite = 0, kBlack = 1, kGray = 920, kRed = 632, kGreen = 416, kBlue = 600,
    kYellow = 400, kMagenta = 616, kCyan = 432, kOrange = 800, kSpring = 820,
    kTeal = 840, kAzure = 860, kViolet = 880, kPink = 900;

var symbols_map = {
   // greek letters
   '#alpha' : '\u03B1',
   '#beta' : '\u03B2',
   '#chi' : '\u03C7',
   '#delta' : '\u03B4',
   '#varepsilon' : '\u03B5',
   '#phi' : '\u03C6',
   '#gamma' : '\u03B3',
   '#eta' : '\u03B7',
   '#iota' : '\u03B9',
   '#varphi' : '\u03C6',
   '#kappa' : '\u03BA',
   '#lambda' : '\u03BB',
   '#mu' : '\u03BC',
   '#nu' : '\u03BD',
   '#omicron' : '\u03BF',
   '#pi' : '\u03C0',
   '#theta' : '\u03B8',
   '#rho' : '\u03C1',
   '#sigma' : '\u03C3',
   '#tau' : '\u03C4',
   '#upsilon' : '\u03C5',
   '#varomega' : '\u03D6',
   '#omega' : '\u03C9',
   '#xi' : '\u03BE',
   '#psi' : '\u03C8',
   '#zeta' : '\u03B6',
   '#Alpha' : '\u0391',
   '#Beta' : '\u0392',
   '#Chi' : '\u03A7',
   '#Delta' : '\u0394',
   '#Epsilon' : '\u0395',
   '#Phi' : '\u03A6',
   '#Gamma' : '\u0393',
   '#Eta' : '\u0397',
   '#Iota' : '\u0399',
   '#vartheta' : '\u03D1',
   '#Kappa' : '\u039A',
   '#Lambda' : '\u039B',
   '#Mu' : '\u039C',
   '#Nu' : '\u039D',
   '#Omicron' : '\u039F',
   '#Pi' : '\u03A0',
   '#Theta' : '\u0398',
   '#Rho' : '\u03A1',
   '#Sigma' : '\u03A3',
   '#Tau' : '\u03A4',
   '#Upsilon' : '\u03A5',
   '#varsigma' : '\u03C2',
   '#Omega' : '\u03A9',
   '#Xi' : '\u039E',
   '#Psi' : '\u03A8',
   '#Zeta' : '\u0396',
   '#varUpsilon' : '\u03D2',
   '#epsilon' : '\u03B5',
   // math symbols

   '#sqrt' : '\u221A',

   // from TLatex tables #2 & #3
   '#leq' : '\u2264',
   '#/' : '\u2044',
   '#infty' : '\u221E',
   '#voidb' : '\u0192',
   '#club' : '\u2663',
   '#diamond' : '\u2666',
   '#heart' : '\u2665',
   '#spade' : '\u2660',
   '#leftrightarrow' : '\u2194',
   '#leftarrow' : '\u2190',
   '#uparrow' : '\u2191',
   '#rightarrow' : '\u2192',
   '#downarrow' : '\u2193',
   '#circ' : '\u02C6', // ^
   '#pm' : '\xB1',
   '#doublequote' : '\u2033',
   '#geq' : '\u2265',
   '#times' : '\xD7',
   '#propto' : '\u221D',
   '#partial' : '\u2202',
   '#bullet' : '\u2022',
   '#divide' : '\xF7',
   '#neq' : '\u2260',
   '#equiv' : '\u2261',
   '#approx' : '\u2248', // should be \u2245 ?
   '#3dots' : '\u2026',
   '#cbar' : '\u007C',
   '#topbar' : '\xAF',
   '#downleftarrow' : '\u21B5',
   '#aleph' : '\u2135',
   '#Jgothic' : '\u2111',
   '#Rgothic' : '\u211C',
   '#voidn' : '\u2118',
   '#otimes' : '\u2297',
   '#oplus' : '\u2295',
   '#oslash' : '\u2205',
   '#cap' : '\u2229',
   '#cup' : '\u222A',
   '#supseteq' : '\u2287',
   '#supset' : '\u2283',
   '#notsubset' : '\u2284',
   '#subseteq' : '\u2286',
   '#subset' : '\u2282',
   '#int' : '\u222B',
   '#in' : '\u2208',
   '#notin' : '\u2209',
   '#angle' : '\u2220',
   '#nabla' : '\u2207',
   '#oright' : '\xAE',
   '#ocopyright' : '\xA9',
   '#trademark' : '\u2122',
   '#prod' : '\u220F',
   '#surd' : '\u221A',
   '#upoint' : '\u22C5',
   '#corner' : '\xAC',
   '#wedge' : '\u2227',
   '#vee' : '\u2228',
   '#Leftrightarrow' : '\u21D4',
   '#Leftarrow' : '\u21D0',
   '#Uparrow' : '\u21D1',
   '#Rightarrow' : '\u21D2',
   '#Downarrow' : '\u21D3',
   '#LT' : '\x3C',
   '#void1' : '\xAE',
   '#copyright' : '\xA9',
   '#void3' : '\u2122',
   '#sum' : '\u2211',
   '#arctop' : '',
   '#lbar' : '',
   '#arcbottom' : '',
   '#void8' : '',
   '#bottombar' : '\u230A',
   '#arcbar' : '',
   '#ltbar' : '',
   '#AA' : '\u212B',
   '#aa' : '\u00E5',
   '#void06' : '',
   '#GT' : '\x3E',
   '#forall' : '\u2200',
   '#exists' : '\u2203',
   '#bar' : '',
   '#vec' : '',
   '#dot' : '\u22C5',
   '#hat' : '\xB7',
   '#ddot' : '',
   '#acute' : '\acute',
   '#grave' : '',
   '#check' : '\u2713',
   '#tilde' : '\u02DC',
   '#slash' : '\u2044',
   '#hbar' : '\u0127',
   '#box' : '',
   '#Box' : '',
   '#parallel' : '',
   '#perp' : '\u22A5',
   '#odot' : ''
};

var tooltip = function() {
   var id = 'tt';
   var top = 3;
   var left = 3;
   var maxw = 150;
   var speed = 10;
   var timer = 20;
   var endalpha = 95;
   var alpha = 0;
   var tt,t,c,b,h;
   var ie = document.all ? true : false;
   return {
      show: function(v, w) {
         if (tt == null) {
            tt = document.createElement('div');
            tt.setAttribute('id',id);
            t = document.createElement('div');
            t.setAttribute('id',id + 'top');
            c = document.createElement('div');
            c.setAttribute('id',id + 'cont');
            b = document.createElement('div');
            b.setAttribute('id',id + 'bot');
            tt.appendChild(t);
            tt.appendChild(c);
            tt.appendChild(b);
            document.body.appendChild(tt);
            tt.style.opacity = 0;
            tt.style.filter = 'alpha(opacity=0)';
            document.onmousemove = this.pos;
         }
         tt.style.display = 'block';
         c.innerHTML = v;
         tt.style.width = w ? w + 'px' : 'auto';
         tt.style.width = 'auto'; // let it be automatically resizing...
         if (!w && ie) {
            t.style.display = 'none';
            b.style.display = 'none';
            tt.style.width = tt.offsetWidth;
            t.style.display = 'block';
            b.style.display = 'block';
         }
         //if (tt.offsetWidth > maxw) { tt.style.width = maxw + 'px'; }
         h = parseInt(tt.offsetHeight) + top;
         clearInterval(tt.timer);
         tt.timer = setInterval( function() { tooltip.fade(1) }, timer );
      },
      pos: function(e) {
         var u = ie ? event.clientY + document.documentElement.scrollTop : e.pageY;
         var l = ie ? event.clientX + document.documentElement.scrollLeft : e.pageX;
         tt.style.top = u + 15 + 'px';//(u - h) + 'px';
         tt.style.left = (l + left) + 'px';
      },
      fade: function(d) {
         var a = alpha;
         if ((a != endalpha && d == 1) || (a != 0 && d == -1)) {
            var i = speed;
            if (endalpha - a < speed && d == 1) {
               i = endalpha - a;
            }
            else if (alpha < speed && d == -1) {
               i = a;
            }
            alpha = a + (i * d);
            tt.style.opacity = alpha * .01;
            tt.style.filter = 'alpha(opacity=' + alpha + ')';
         }
         else {
            clearInterval(tt.timer);
            if (d == -1) { tt.style.display = 'none'; }
         }
      },
      hide: function() {
         if (tt == null) return;
         clearInterval(tt.timer);
         tt.timer = setInterval( function() { tooltip.fade(-1) }, timer );
      }
   };
}();

/**
 * @author alteredq / http://alteredqualia.com/
 * @author mr.doob / http://mrdoob.com/
 */
var Detector = {
   canvas: !! window.CanvasRenderingContext2D,
   webgl: ( function () { try { return !! window.WebGLRenderingContext && !! document.createElement( 'canvas' ).getContext( 'experimental-webgl' ); } catch( e ) { return false; } } )(),
   workers: !! window.Worker, fileapi: window.File && window.FileReader && window.FileList && window.Blob
};

/*
 * Function that generates all root colors
 */
function generateAllColors () {
   var colorMap = new Array(
      'rgb(255, 255, 255)',
      'rgb(0, 0, 0)',
      'rgb(255, 0, 0)',
      'rgb(0, 255, 0)',
      'rgb(0, 0, 255)',
      'rgb(255, 255, 0)',
      'rgb(255, 0, 255)',
      'rgb(0, 255, 255)',
      'rgb(89, 211, 84)',
      'rgb(89, 84, 216)',
      'rgb(254, 254, 254)',
      'rgb(191, 181, 173)',
      'rgb(76, 76, 76)',
      'rgb(102, 102, 102)',
      'rgb(127, 127, 127)',
      'rgb(153, 153, 153)',
      'rgb(178, 178, 178)',
      'rgb(204, 204, 204)',
      'rgb(229, 229, 229)',
      'rgb(242, 242, 242)',
      'rgb(204, 198, 170)',
      'rgb(204, 198, 170)',
      'rgb(193, 191, 168)',
      'rgb(186, 181, 163)',
      'rgb(178, 165, 150)',
      'rgb(183, 163, 155)',
      'rgb(173, 153, 140)',
      'rgb(155, 142, 130)',
      'rgb(135, 102, 86)',
      'rgb(175, 206, 198)',
      'rgb(132, 193, 163)',
      'rgb(137, 168, 160)',
      'rgb(130, 158, 140)',
      'rgb(173, 188, 198)',
      'rgb(122, 142, 153)',
      'rgb(117, 137, 145)',
      'rgb(104, 130, 150)',
      'rgb(109, 122, 132)',
      'rgb(124, 153, 209)',
      'rgb(127, 127, 155)',
      'rgb(170, 165, 191)',
      'rgb(211, 206, 135)',
      'rgb(221, 186, 135)',
      'rgb(188, 158, 130)',
      'rgb(198, 153, 124)',
      'rgb(191, 130, 119)',
      'rgb(206, 94, 96)',
      'rgb(170, 142, 147)',
      'rgb(165, 119, 122)',
      'rgb(147, 104, 112)',
      'rgb(211, 89, 84)');

   var circleColors = [632, 416, 600, 400, 616, 432];

   var rectangleColors = [800, 820, 840, 860, 880, 900];

   var set1 = [ 255,204,204, 255,153,153, 204,153,153, 255,102,102, 204,102,102
               ,153,102,102, 255, 51, 51, 204, 51, 51, 153, 51, 51, 102, 51, 51
               ,255,  0,  0, 204,  0,  0, 153,  0,  0, 102,  0,  0,  51,  0,  0];
   var set2 = [ 204,255,204, 153,255,153, 153,204,153, 102,255,102, 102,204,102
               ,102,153,102, 51,255, 51,  51,204, 51,  51,153, 51,  51,102, 51
               ,  0,255,  0,   0,204,  0,   0,153,  0,   0,102,  0,  0, 51,  0];
   var set3 = [ 204,204,255, 153,153,255, 153,153,204, 102,102,255, 102,102,204
               ,102,102,153,  51, 51,255,  51, 51,204,  51, 51,153,  51, 51,102
               ,  0,  0,255,   0,  0,204,   0,  0,153,   0,  0,102,   0,  0, 51];
   var set4 = [ 255,255,204, 255,255,153, 204,204,153, 255,255,102, 204,204,102
               ,153,153,102, 255,255, 51, 204,204, 51, 153,153, 51, 102,102, 51
               ,255,255,  0, 204,204,  0, 153,153,  0, 102,102,  0,  51, 51,  0];
   var set5 = [ 255,204,255, 255,153,255, 204,153,204, 255,102,255, 204,102,204
               ,153,102,153, 255, 51,255, 204, 51,204, 153, 51,153, 102, 51,102
               ,255,  0,255, 204,  0,204, 153,  0,153, 102,  0,102,  51,  0, 51];
   var set6 = [ 204,255,255, 153,255,255, 153,204,204, 102,255,255, 102,204,204
               ,102,153,153,  51,255,255,  51,204,204,  51,153,153,  51,102,102
               ,  0,255,255,   0,204,204,   0,153,153,   0,102,102,   0, 51,  51];

   var circleSets = new Array(set1, set2, set3, set4, set5, set6);

   var set7 = [ 255,204,153,  204,153,102,  153,102, 51,  153,102,  0,  204,153, 51
                ,255,204,102,  255,153,  0,  255,204, 51,  204,153,  0,  255,204,  0
                ,255,153, 51,  204,102,  0,  102, 51,  0,  153, 51,  0,  204,102, 51
                ,255,153,102,  255,102,  0,  255,102, 51,  204, 51,  0,  255, 51,  0];
   var set8 = [ 153,255, 51,  102,204,  0,   51,102,  0,   51,153,  0,  102,204, 51
               ,153,255,102,  102,255,  0,  102,255, 51,   51,204,  0,   51,255, 0
               ,204,255,153,  153,204,102,  102,153, 51,  102,153,  0,  153,204, 51
               ,204,255,102,  153,255,  0,  204,255, 51,  153,204,  0,  204,255,  0];
   var set9 = [ 153,255,204,  102,204,153,   51,153,102,    0,153,102,   51,204,153
               ,102,255,204,    0,255,102,   51,255,204,    0,204,153,    0,255,204
               , 51,255,153,    0,204,102,    0,102, 51,    0,153, 51,   51,204,102
               ,102,255,153,    0,255,153,   51,255,102,    0,204, 51,    0,255, 51];
   var set10 = [153,204,255,  102,153,204,   51,102,153,    0, 51,153,   51,102,204
               ,102,153,255,    0,102,255,   51,102,255,    0, 51,204,    0, 51,255
               , 51,153,255,    0,102,204,    0, 51,102,    0,102,153,   51,153,204
               ,102,204,255,    0,153,255,   51,204,255,    0,153,204,    0,204,255];
   var set11 = [204,153,255,  153,102,204,  102, 51,153,  102,  0,153,  153, 51,204
               ,204,102,255,  153,  0,255,  204, 51,255,  153,  0,204,  204,  0,255
               ,153, 51,255,  102,  0,204,   51,  0,102,   51,  0,153,  102, 51,204
               ,153,102,255,  102,  0,255,  102, 51,255,   51,  0,204,   51,  0,255];
   var set12 = [255, 51,153,  204,  0,102,  102,  0, 51,  153,  0, 51,  204, 51,102
               ,255,102,153,  255,  0,102,  255, 51,102,  204,  0, 51,  255,  0, 51
               ,255,153,204,  204,102,153,  153, 51,102,  153,  0,102,  204, 51,153
               ,255,102,204,  255,  0,153,  204,  0,153,  255, 51,204,  255,  0,153];

   var rectSets = new Array(set7, set8, set9, set10, set11, set12);

   /*
    * Define circle colors
    */
   for(var i = 0; i < 6; i++) {
      for(var j = 0; j < 15; j++) {
         var colorn = circleColors[i] + j - 10;
         colorMap[colorn] = 'rgb(' + circleSets[i][3*j] + ', ' + circleSets[i][3*j+1] + ', ' + circleSets[i][3*j+2] + ')';
         colorn = rectangleColors[i] + j - 9;
         colorMap[colorn] = 'rgb('+ rectSets[i][3*j] + ', ' + rectSets[i][3*j+1] + ', ' + rectSets[i][3*j+2] + ')';
      }
    }
    return colorMap;
};

function getFontDetails(fontName) {
   var weight = "";
   var style = "";
   var name = "Arial";

   if (fontName.indexOf("bold") != -1) {
      weight = "bold";
      //The first 5 characters are removed because "bold " is always first when it occurs
      fontName = fontName.substring(5, fontName.length);
   }
   if (fontName.charAt(0) == 'i') {
      style = "italic";
      fontName = fontName.substring(7, fontName.length);
   }
   else if (fontName.charAt(0) == 'o') {
      style = "oblique";
      fontName = fontName.substring(8, fontName.length);
   }
   if (name == 'Symbol') {
      weight = "";
      style = "";
   }
   return {
      'weight' : weight,
      'style'  : style,
      'name'   : fontName
   };
};

/*
 * Function that returns the SVG symbol type identifier for a given root matker
 * The result is an array with 3 elements:
 *    the first is the identifier of the root marker in the SVG symbols
 *    the second is true if the shape is filled and false if it is open
 *    the third is true if the shape should be rotated
 * The identifier will be 6 if the shape is a star or 7 if it is '*'
 */
function getRootMarker(markers, i) {
   var marker = markers[i];
   var shape = 0;
   var toFill = true;
   var toRotate = false;

   if (typeof(marker) != 'undefined') {
      var fst = marker.charAt(0);
      switch (fst) {
         case 'd':
            shape = 7;
            return {'shape' : shape};
         case 'o':
            toFill = false;
            break;
         case 'g':
            toRotate = true;
      }

      var type = marker.substr(1, marker.length);
      switch (type) {
         case "circle":
            shape = 0;
            break;
         case "cross":
            shape = 1;
            break;
         case "diamond":
            shape = 2;
            break;
         case "square":
            shape = 3;
            break;
         case "triangle-up":
            shape = 4;
            break;
         case "triangle-down":
            shape = 5;
            break;
         case "star":
            shape = 6;
            break;
      };
   }
   return {
      'shape'    : shape,
      'toFill'   : toFill,
      'toRotate' : toRotate
   };
};

function stringWidth(svg, line, font_name, font_weight, font_size, font_style) {
   /* compute the bounding box of a string by using temporary svg:text */
   var text = svg.append("svg:text")
       .attr("class", "temp_text")
       .attr("font-family", font_name)
       .attr("font-weight", font_weight)
       .attr("font-style", font_style)
       .attr("font-size", font_size)
       .style("opacity", 0)
       .text(line);
   w = text.node().getBBox().width;
   text.remove();
   return w;
}

function format_id(id) {
   /* format the string id to remove specials characters
      (that cannot be used in id strings) */
   var g_id = id;
   if (g_id == "") g_id = "random_histo_" + random_id++;
   while (g_id.indexOf(' ') != -1)
      g_id = g_id.replace(' ', '_');
   while (g_id.indexOf(':') != -1)
      g_id = g_id.replace(':', '_');
   while (g_id.indexOf('.') != -1)
      g_id = g_id.replace('.', '_');
   while (g_id.indexOf('>') != -1)
      g_id = g_id.replace('>', 'gt');
   while (g_id.indexOf('<') != -1)
      g_id = g_id.replace('<', 'lt');
   while (g_id.indexOf('\\') != -1)
      g_id = g_id.replace('\\', '');
   while (g_id.indexOf('\'') != -1)
      g_id = g_id.replace('\'', '');
   while (g_id.indexOf('(') != -1)
      g_id = g_id.replace('(', '_');
   while (g_id.indexOf(')') != -1)
      g_id = g_id.replace(')', '_');
   while (g_id.indexOf('/') != -1)
      g_id = g_id.replace('/', '_');
   while (g_id.indexOf('-') != -1)
      g_id = g_id.replace('-', '_');
   while (g_id.indexOf('[') != -1)
      g_id = g_id.replace('[', '_');
   while (g_id.indexOf(']') != -1)
      g_id = g_id.replace(']', '_');
   return g_id;
};

/*
    Polyfill for touch dblclick
    http://mckamey.mit-license.org
*/
function doubleTap(elem, speed, distance) {
   if (!('ontouchstart' in elem)) {
      // non-touch has native dblclick and no need for polyfill
      return;
   }
   // default dblclick speed to half sec
   speed = Math.abs(+speed) || 500;//ms
   // default dblclick distance to within 40x40 area
   distance = Math.abs(+distance) || 40;//px

   var taps, x, y;
   var reset = function() {
      // reset state
      taps = 0;
      x = NaN;
      y = NaN;
   };
   reset();

   elem.addEventListener('touchstart', function(e) {
      var touch = e.changedTouches[0] || {}, oldX = x, oldY = y;

      taps++;
      x = +touch.pageX || +touch.clientX || +touch.screenX;
      y = +touch.pageY || +touch.clientY || +touch.screenY;

      // NaN will always be false
      if (Math.abs(oldX-x) < distance && Math.abs(oldY-y) < distance) {
         // fire dblclick event
         var e2 = document.createEvent('MouseEvents');
         if (e2.initMouseEvent) {
             e2.initMouseEvent(
                 'dblclick',
                 true,                   // dblclick bubbles
                 true,                   // dblclick cancelable
                 e.view,                 // copy view
                 taps,                   // click count
                 touch.screenX,          // copy coordinates
                 touch.screenY,
                 touch.clientX,
                 touch.clientY,
                 e.ctrlKey,              // copy key modifiers
                 e.altKey,
                 e.shiftKey,
                 e.metaKey,
                 e.button,               // copy button 0: left, 1: middle, 2: right
                 touch.target);          // copy target
         }
         elem.dispatchEvent(e2);
      }
      setTimeout(reset, speed);
   }, false);

   elem.addEventListener('touchmove', function(e) {
      reset();
   }, false);
};

function createFillPatterns(svg, id, color) {
   // create fill patterns - only if they don't exists yet
   var line_color = JSROOTPainter.getRootColor(color);
   for (var i=0; i<defs[0][0]['childNodes'].length;++i) {
      if (defs[0][0]['childNodes'][i]['id'] == "pat" + id + "_" + color &&
          defs[0][0]['childNodes'][i]['style']['stroke'] == line_color)
         return;
   }
   switch (id) {
      case 3001:
         defs.append('svg:pattern')
            .attr("id", "pat3001_"+color)
            .attr("patternUnits", "userSpaceOnUse")
            .attr("width", "3px")
            .attr("height", "2px")
            .style("stroke", line_color)
            .append('svg:rect')
            .attr("x", 0)
            .attr("y", 0)
            .attr("width", 1)
            .attr("height", 1)
            .style("stroke", line_color)
            .append('svg:rect')
            .attr("x", 2)
            .attr("y", 0)
            .attr("width", 1)
            .attr("height", 1)
            .style("stroke", line_color)
            .append('svg:rect')
            .attr("x", 1)
            .attr("y", 1)
            .attr("width", 1)
            .attr("height", 1)
            .style("stroke", line_color);
         break;
      case 3002:
         defs.append('svg:pattern')
            .attr("id", "pat3002_"+color)
            .attr("patternUnits", "userSpaceOnUse")
            .attr("width", "4px")
            .attr("height", "2px")
            .style("stroke", line_color)
            .append('svg:rect')
            .attr("x", 1)
            .attr("y", 0)
            .attr("width", 1)
            .attr("height", 1)
            .style("stroke", line_color)
            .append('svg:rect')
            .attr("x", 3)
            .attr("y", 1)
            .attr("width", 1)
            .attr("height", 1)
            .style("stroke", line_color);
         break;
      case 3003:
         defs.append('svg:pattern')
            .attr("id", "pat3003_"+color)
            .attr("patternUnits", "userSpaceOnUse")
            .attr("width", "4px")
            .attr("height", "4px")
            .style("stroke", line_color)
            .append('svg:rect')
            .attr("x", 2)
            .attr("y", 1)
            .attr("width", 1)
            .attr("height", 1)
            .style("stroke", line_color)
            .append('svg:rect')
            .attr("x", 0)
            .attr("y", 3)
            .attr("width", 1)
            .attr("height", 1)
            .style("stroke", line_color);
         break;
      case 3004:
         defs.append('svg:pattern')
            .attr("id", "pat3004_"+color)
            .attr("patternUnits", "userSpaceOnUse")
            .attr("width", "8px")
            .attr("height", "8px")
            .style("stroke", line_color)
            .append("svg:line")
            .attr("x1", 8)
            .attr("y1", 0)
            .attr("x2", 0)
            .attr("y2", 8)
            .style("stroke", line_color)
            .style("stroke-width", 1);
         break;
      case 3005:
         defs.append('svg:pattern')
            .attr("id", "pat3005_"+color)
            .attr("patternUnits", "userSpaceOnUse")
            .attr("width", "8px")
            .attr("height", "8px")
            .style("stroke", line_color)
            .append("svg:line")
            .attr("x1", 0)
            .attr("y1", 0)
            .attr("x2", 8)
            .attr("y2", 8)
            .style("stroke", line_color)
            .style("stroke-width", 1);
         break;
      case 3006:
         defs.append('svg:pattern')
            .attr("id", "pat3006_"+color)
            .attr("patternUnits", "userSpaceOnUse")
            .attr("width", "4px")
            .attr("height", "4px")
            .style("stroke", line_color)
            .append("svg:line")
            .attr("x1", 1)
            .attr("y1", 0)
            .attr("x2", 1)
            .attr("y2", 3)
            .style("stroke", line_color)
            .style("stroke-width", 1);
         break;
      case 3007:
         defs.append('svg:pattern')
            .attr("id", "pat3007_"+color)
            .attr("patternUnits", "userSpaceOnUse")
            .attr("width", "4px")
            .attr("height", "4px")
            .style("stroke", line_color)
            .append("svg:line")
            .attr("x1", 0)
            .attr("y1", 1)
            .attr("x2", 3)
            .attr("y2", 1)
            .style("stroke", line_color)
            .style("stroke-width", 1);
         break;
      default: /* == 3004 */
         defs.append('svg:pattern')
            .attr("id", "pat"+id+"_"+color)
            .attr("patternUnits", "userSpaceOnUse")
            .attr("width", "8px")
            .attr("height", "8px")
            .style("stroke", line_color)
            .append("svg:line")
            .attr("x1", 8)
            .attr("y1", 0)
            .attr("x2", 0)
            .attr("y2", 8)
            .style("stroke", line_color)
            .style("stroke-width", 1);
         break;
   }
};

(function(){

   if (typeof JSROOTPainter == 'object'){
      var e1 = new Error('JSROOTPainter is already defined');
      e1.source = 'JSROOTPainter.js';
      throw e1;
   }

   if (typeof d3 != 'object') {
      var e1 = new Error('This extension requires d3.js.js');
      e1.source = 'JSROOTPainter.js';
      throw e1;
   }

   // Initialize Custom colors
   var root_colors = generateAllColors();

   // Initialize colors of the default palette
   var default_palette = new Array();

   //Initialize ROOT markers
   var root_markers = new Array('fcircle','fcircle', 'fcross', 'dcross', 'ocircle',
      'gcross', 'fcircle', 'fcircle', 'fcircle', 'fcircle', 'fcircle',
      'fcircle', 'fcircle', 'fcircle', 'fcircle', 'fcircle', 'fcircle',
      'fcircle', 'fcircle', 'fcircle', 'fcircle', 'fsquare', 'ftriangle-up',
      'ftriangle-down', 'ocircle', 'osquare', 'otriangle-up', 'odiamond',
      'ocross', 'fstar', 'ostar', 'dcross', 'otriangle-down', 'fdiamond',
      'fcross');

   var root_fonts = new Array('Arial', 'Times New Roman',
      'bold Times New Roman', 'bold italic Times New Roman',
      'Arial', 'oblique Arial', 'bold Arial', 'bold oblique Arial',
      'Courier New', 'oblique Courier New', 'bold Courier New',
      'bold oblique Courier New', 'Symbol', 'Times New Roman',
      'Wingdings', 'Symbol');

   var root_line_styles = new Array("", "", "3, 3", "1, 2", "3, 4, 1, 4",
         "5, 3, 1, 3", "5, 3, 1, 3, 1, 3, 1, 3", "5, 5",
         "5, 3, 1, 3, 1, 3", "20, 5", "20, 10, 1, 10", "1, 2");

   JSROOTPainter = {};

   JSROOTPainter.version = '2.0 2012/11/10';

   /*
    * Helper functions
    */

   JSROOTPainter.clearCuts = function(chopt) {
      /* decode string "chopt" and remove graphical cuts */
      var left = chopt.indexOf('[');
      if (left == -1) return chopt;
      var right = chopt.indexOf(']');
      if (right == -1) return chopt;
      var nch = right-left;
      if (nch < 2) return chopt;
      for (i=0;i<=nch;i++) chopt[left+i] = ' ';
      return chopt;
   };

   JSROOTPainter.decodeOptions = function(opt, histo, pad) {
      /* decode string 'opt' and fill the option structure */
      var hdim = 1; // histo['fDimension'];
      if (histo['_typename'].match(/\bJSROOTIO.TH2/)) hdim = 2;
      if (histo['_typename'].match(/\bJSROOTIO.TH3/)) hdim = 3;
      var nch = opt.length;
      var option = { 'Axis': 0, 'Bar': 0, 'Curve': 0, 'Error': 0, 'Hist': 0,
         'Line': 0, 'Mark': 0, 'Fill': 0, 'Same': 0, 'Func': 0, 'Scat': 0,
         'Star': 0, 'Arrow': 0, 'Box': 0, 'Text': 0, 'Char': 0, 'Color': 0,
         'Contour': 0, 'Logx': 0, 'Logy': 0, 'Logz': 0, 'Lego': 0, 'Surf': 0,
         'Off': 0, 'Tri': 0, 'Proj': 0, 'AxisPos': 0, 'Spec': 0, 'Pie': 0,
         'List': 0, 'Zscale': 0, 'FrontBox': 1, 'BackBox': 1, 'System': kCARTESIAN,
         'HighRes': 0, 'Zero': 0
      };
      //check for graphical cuts
      var chopt = opt.toUpperCase();
      chopt = JSROOTPainter.clearCuts(chopt);
      if (hdim > 1) option.Scat = 1;
      if (!nch) option.Hist = 1;
      if ('fFunctions' in histo && histo['fFunctions'].length > 0) option.Func = 2;
      if ('fSumw2' in histo && histo['fSumw2'].length > 0 && hdim == 1) option.Error = 2;
      var l = chopt.indexOf('SPEC');
      if (l != -1) {
         option.Scat = 0;
         chopt = chopt.replace('SPEC', '    ');
         var bs = 0;
         l = chopt.indexOf('BF(');
         if (l != -1) {
            bs = parseInt(chopt)
         }
         option.Spec = Math.max(1600, bs);
         return option;
      }
      l = chopt.indexOf('GL');
      if (l != -1) {
         chopt = chopt.replace('GL', '  ');
      }
      l = chopt.indexOf('X+');
      if (l != -1) {
         option.AxisPos = 10;
         chopt = chopt.replace('X+', '  ');
      }
      l = chopt.indexOf('Y+');
      if (l != -1) {
         option.AxisPos += 1;
         chopt = chopt.replace('Y+', '  ');
      }
      if ((option.AxisPos == 10 || option.AxisPos == 1) && (nch == 2)) option.Hist = 1;
      if (option.AxisPos == 11 && nch == 4) option.Hist = 1;
      l = chopt.indexOf('SAMES');
      if (l != -1) {
         if (nch == 5) option.Hist = 1;
         option.Same = 2;
         chopt = chopt.replace('SAMES', '     ');
      }
      l = chopt.indexOf('SAME');
      if (l != -1) {
         if (nch == 4) option.Hist = 1;
         option.Same = 1;
         chopt = chopt.replace('SAME', '    ');
      }
      l = chopt.indexOf('PIE');
      if (l != -1) {
         option.Pie = 1;
         chopt = chopt.replace('PIE', '   ');
      }
      l = chopt.indexOf('LEGO');
      if (l != -1) {
         option.Scat = 0;
         option.Lego = 1;
         chopt = chopt.replace('LEGO', '    ');
         if (chopt[l+4] == '1') { option.Lego = 11; chopt[l+4] = ' '; }
         if (chopt[l+4] == '2') { option.Lego = 12; chopt[l+4] = ' '; }
         if (chopt[l+4] == '3') { option.Lego = 13; chopt[l+4] = ' '; }
         l = chopt.indexOf('FB'); if (l != -1) { option.FrontBox = 0; chopt = chopt.replace('FB', '  '); }
         l = chopt.indexOf('BB'); if (l != -1) { option.BackBox = 0;  chopt = chopt.replace('BB', '  '); }
         l = chopt.indexOf('0'); if (l != -1) { option.Zero = 1;  chopt = chopt.replace('0', ' '); }
      }
      l = chopt.indexOf('SURF');
      if (l != -1) {
         option.Scat = 0;
         option.Surf = 1; chopt = chopt.replace('SURF', '    ');
         if (chopt[l+4] == '1') { option.Surf = 11; chopt[l+4] = ' '; }
         if (chopt[l+4] == '2') { option.Surf = 12; chopt[l+4] = ' '; }
         if (chopt[l+4] == '3') { option.Surf = 13; chopt[l+4] = ' '; }
         if (chopt[l+4] == '4') { option.Surf = 14; chopt[l+4] = ' '; }
         if (chopt[l+4] == '5') { option.Surf = 15; chopt[l+4] = ' '; }
         if (chopt[l+4] == '6') { option.Surf = 16; chopt[l+4] = ' '; }
         if (chopt[l+4] == '7') { option.Surf = 17; chopt[l+4] = ' '; }
         l = chopt.indexOf('FB'); if (l != -1) { option.FrontBox = 0; chopt = chopt.replace('FB', '  '); }
         l = chopt.indexOf('BB'); if (l != -1) { option.BackBox = 0; chopt = chopt.replace('BB', '  '); }
      }
      l = chopt.indexOf('TF3');
      if (l != -1) {
         l = chopt.indexOf('FB'); if (l != -1) { option.FrontBox = 0; chopt = chopt.replace('FB', '  '); }
         l = chopt.indexOf('BB'); if (l != -1) { option.BackBox = 0; chopt = chopt.replace('BB', '  '); }
      }
      l = chopt.indexOf('ISO');
      if (l != -1) {
         l = chopt.indexOf('FB'); if (l != -1) { option.FrontBox = 0; chopt = chopt.replace('FB', '  '); }
         l = chopt.indexOf('BB'); if (l != -1) { option.BackBox = 0; chopt = chopt.replace('BB', '  '); }
      }
      l = chopt.indexOf('LIST'); if (l != -1) { option.List = 1; chopt = chopt.replace('LIST', '  '); }
      l = chopt.indexOf('CONT');
      if (l != -1) {
         chopt = chopt.replace('CONT', '    ');
         if (hdim > 1) {
            option.Scat = 0;
            option.Contour = 1;
            if (chopt[l+4] == '1') { option.Contour = 11; chopt[l+4] = ' '; }
            if (chopt[l+4] == '2') { option.Contour = 12; chopt[l+4] = ' '; }
            if (chopt[l+4] == '3') { option.Contour = 13; chopt[l+4] = ' '; }
            if (chopt[l+4] == '4') { option.Contour = 14; chopt[l+4] = ' '; }
            if (chopt[l+4] == '5') { option.Contour = 15; chopt[l+4] = ' '; }
         } else {
            option.Hist = 1;
         }
      }
      l = chopt.indexOf('HBAR');
      if (l != -1) {
         option.Hist = 0;
         option.Bar = 20; chopt = chopt.replace('HBAR', '    ');
         if (chopt[l+4] == '1') { option.Bar = 21; chopt[l+4] = ' '; }
         if (chopt[l+4] == '2') { option.Bar = 22; chopt[l+4] = ' '; }
         if (chopt[l+4] == '3') { option.Bar = 23; chopt[l+4] = ' '; }
         if (chopt[l+4] == '4') { option.Bar = 24; chopt[l+4] = ' '; }
      }
      l = chopt.indexOf('BAR');
      if (l != -1) {
         option.Hist = 0;
         option.Bar = 10; chopt = chopt.replace('BAR', '   ');
         if (chopt[l+3] == '1') { option.Bar = 11; chopt[l+3] = ' '; }
         if (chopt[l+3] == '2') { option.Bar = 12; chopt[l+3] = ' '; }
         if (chopt[l+3] == '3') { option.Bar = 13; chopt[l+3] = ' '; }
         if (chopt[l+3] == '4') { option.Bar = 14; chopt[l+3] = ' '; }
      }
      l = chopt.indexOf('ARR' );
      if (l != -1) {
         chopt = chopt.replace('ARR', '   ');
         if (hdim > 1) {
            option.Arrow  = 1;
            option.Scat = 0;
         } else {
            option.Hist = 1;
         }
      }
      l = chopt.indexOf('BOX' );
      if (l != -1) {
         chopt = chopt.replace('BOX', '   ');
         if (hdim>1) {
            Hoption.Scat = 0;
            Hoption.Box  = 1;
            if (chopt[l+3] == '1') { option.Box = 11; chopt[l+3] = ' '; }
         } else {
            option.Hist = 1;
         }
      }
      l = chopt.indexOf('COLZ');
      if (l != -1) {
         chopt = chopt.replace('COLZ', '');
         if (hdim>1) {
            option.Color  = 2;
            option.Scat   = 0;
            option.Zscale = 1;
         } else {
            option.Hist = 1;
         }
      }
      l = chopt.indexOf('COL' );
      if (l != -1) {
         chopt = chopt.replace('COL', '   ');
         if (hdim>1) {
            option.Color = 1;
            option.Scat  = 0;
         } else {
            option.Hist = 1;
         }
      }
      l = chopt.indexOf('CHAR'); if (l != -1) { option.Char = 1; chopt = chopt.replace('CHAR', '    '); option.Scat = 0; }
      l = chopt.indexOf('FUNC'); if (l != -1) { option.Func = 2; chopt = chopt.replace('FUNC', '    '); option.Hist = 0; }
      l = chopt.indexOf('HIST'); if (l != -1) { option.Hist = 2; chopt = chopt.replace('HIST', '    '); option.Func = 0; option.Error = 0; }
      l = chopt.indexOf('AXIS'); if (l != -1) { option.Axis = 1; chopt = chopt.replace('AXIS', '    '); }
      l = chopt.indexOf('AXIG'); if (l != -1) { option.Axis = 2; chopt = chopt.replace('AXIG', '    '); }
      l = chopt.indexOf('SCAT'); if (l != -1) { option.Scat = 1; chopt = chopt.replace('SCAT', '    '); }
      l = chopt.indexOf('TEXT');
      if (l != -1) {
         var angle = parseInt(chopt);
         if (!isNaN(angle)) {
            if (angle < 0)  angle = 0;
            if (angle > 90) angle = 90;
            option.Text = 1000 + angle;
         }
         else {
            option.Text = 1;
         }
         chopt = chopt.replace('TEXT', '    ');
         l = chopt.indexOf('N');
         if (l != -1 && histo['_typename'].match(/\bJSROOTIO.TH2Poly/)) option.Text += 3000;
         option.Scat = 0;
      }
      l = chopt.indexOf('POL');  if (l != -1) { option.System = kPOLAR;       chopt = chopt.replace('POL', '   '); }
      l = chopt.indexOf('CYL');  if (l != -1) { option.System = kCYLINDRICAL; chopt = chopt.replace('CYL', '   '); }
      l = chopt.indexOf('SPH');  if (l != -1) { option.System = kSPHERICAL;   chopt = chopt.replace('SPH', '   '); }
      l = chopt.indexOf('PSR');  if (l != -1) { option.System = kRAPIDITY;    chopt = chopt.replace('PSR', '   '); }
      l = chopt.indexOf('TRI');
      if (l != -1) {
         option.Scat = 0;
         option.Color  = 0;
         option.Tri = 1; chopt = chopt.replace('TRI', '   ');
         l = chopt.indexOf('FB');   if (l != -1) { option.FrontBox = 0; chopt = chopt.replace('FB', '  '); }
         l = chopt.indexOf('BB');   if (l != -1) { option.BackBox = 0;  chopt = chopt.replace('BB', '  '); }
         l = chopt.indexOf('ERR');  if (l != -1) chopt = chopt.replace('ERR', '   ');
      }
      l = chopt.indexOf('AITOFF');
      if (l != -1) {
         Hoption.Proj = 1; chopt = chopt.replace('AITOFF', '      ');  // Aitoff projection
      }
      l = chopt.indexOf('MERCATOR');
      if (l != -1) {
         option.Proj = 2; chopt = chopt.replace('MERCATOR', '        ');  // Mercator projection
      }
      l = chopt.indexOf('SINUSOIDAL');
      if (l != -1) {
         option.Proj = 3; chopt = chopt.replace('SINUSOIDAL', '         ');  // Sinusoidal projection
      }
      l = chopt.indexOf('PARABOLIC');
      if (l != -1) {
         option.Proj = 4; chopt = chopt.replace('PARABOLIC', '         ');  // Parabolic projection
      }
      if (option.Proj > 0) {
         option.Scat = 0;
         option.Contour = 14;
      }
      if (chopt.indexOf('A') != -1)    option.Axis = -1;
      if (chopt.indexOf('B') != -1)    option.Bar  = 1;
      if (chopt.indexOf('C') != -1)  { option.Curve =1; option.Hist = -1; }
      if (chopt.indexOf('F') != -1)    option.Fill =1;
      if (chopt.indexOf('][') != -1) { option.Off  =1; option.Hist =1; }
      if (chopt.indexOf('F2') != -1)   option.Fill =2;
      if (chopt.indexOf('L') != -1)  { option.Line =1; option.Hist = -1; }
      if (chopt.indexOf('P') != -1)  { option.Mark =1; option.Hist = -1; }
      if (chopt.indexOf('Z') != -1)    option.Zscale =1;
      if (chopt.indexOf('*') != -1)    option.Star =1;
      if (chopt.indexOf('H') != -1)    option.Hist =2;
      if (chopt.indexOf('P0') != -1)   option.Mark =10;
      if (histo['_typename'].match(/\bJSROOTIO.TH2Poly/)) {
         if (option.Fill + option.Line + option.Mark != 0 ) option.Scat = 0;
      }
      if (chopt.indexOf('E') != -1) {
         if (hdim == 1) {
            option.Error = 1;
            if (chopt.indexOf('E0') != -1)  option.Error = 10;
            if (chopt.indexOf('E1') != -1)  option.Error = 11;
            if (chopt.indexOf('E2') != -1)  option.Error = 12;
            if (chopt.indexOf('E3') != -1)  option.Error = 13;
            if (chopt.indexOf('E4') != -1)  option.Error = 14;
            if (chopt.indexOf('E5') != -1)  option.Error = 15;
            if (chopt.indexOf('E6') != -1)  option.Error = 16;
            if (chopt.indexOf('X0') != -1) {
               if (option.Error == 1)  option.Error += 20;
               option.Error += 10;
            }
            if (option.Text && histo['_typename'].match(/\bJSROOTIO.TProfile/)) {
               option.Text += 2000;
               option.Error = 0;
            }
         } else {
            if (option.Error == 0) {
               option.Error = 100;
               option.Scat  = 0;
            }
            if (option.Text) {
               option.Text += 2000;
               option.Error = 0;
            }
         }
      }
      if (chopt.indexOf('9') != -1)  option.HighRes = 1;
      if (option.Surf == 15) {
         if (option.System == kPOLAR || option.System == kCARTESIAN) {
            option.Surf = 13;
            //Warning('MakeChopt','option SURF5 is not supported in Cartesian and Polar modes');
         }
      }
      if (pad && typeof(pad) != 'undefined') {
         // Copy options from current pad
         option.Logx = pad['fLogx'];
         option.Logy = pad['fLogy'];
         option.Logz = pad['fLogz'];
      }
      //  Check options incompatibilities
      if (option.Bar == 1) option.Hist = -1;
      return option;
   };

   JSROOTPainter.getRootColor = function(color) {
      return root_colors[color];
   };

   JSROOTPainter.padtoX = function(pad, x) {
      // Convert x from pad to X.
      if (pad['fLogx'] && x < 50) return Math.exp(2.302585092994 * x);
      return x;
   };

   JSROOTPainter.padtoY = function(pad, y) {
      // Convert y from pad to Y.
      if (pad['fLogy'] && y < 50) return Math.exp(2.302585092994 * y);
      return y;
   };

   JSROOTPainter.xtoPad = function(x, pad) {
      if (pad['fLogx']) {
         if (x > 0)
            x = JSROOTMath.log10(x);
         else
            x = pad['fUxmin'];
      }
      return x;
   };

   JSROOTPainter.ytoPad = function(y, pad) {
      if (pad['fLogy']) {
         if (y > 0)
            y = JSROOTMath.log10(y);
         else
            y = pad['fUymin'];
      }
      return y;
   };

   /**
    * Converts an HSL color value to RGB. Conversion formula
    * adapted from http://en.wikipedia.org/wiki/HSL_color_space.
    * Assumes h, s, and l are contained in the set [0, 1] and
    * returns r, g, and b in the set [0, 255].
    *
    * @param   Number  h       The hue
    * @param   Number  s       The saturation
    * @param   Number  l       The lightness
    * @return  Array           The RGB representation
    */
   JSROOTPainter.HLStoRGB = function(h, l, s) {
      var r, g, b;
      if (s < 1e-300) {
         r = g = b = l; // achromatic
      } else {
         function hue2rgb(p, q, t){
            if (t < 0) t += 1;
            if (t > 1) t -= 1;
            if (t < 1/6) return p + (q - p) * 6 * t;
            if (t < 1/2) return q;
            if (t < 2/3) return p + (q - p) * (2/3 - t) * 6;
            return p;
         }
         var q = l < 0.5 ? l * (1 + s) : l + s - l * s;
         var p = 2 * l - q;
         r = hue2rgb(p, q, h + 1/3);
         g = hue2rgb(p, q, h);
         b = hue2rgb(p, q, h - 1/3);
      }
      return 'rgb('+Math.round(r * 255)+', '+Math.round(g * 255)+', '+Math.round(b * 255)+')';
   };

   JSROOTPainter.getMinMax = function(hist, what) {
      if (what == 'max' && hist['fMaximum'] != -1111) return hist['fMaximum'];
      if (what == 'min' && hist['fMinimum'] != -1111) return hist['fMinimum'];
      var bin, binx, biny, binz;
      var xfirst  = 1;;
      var xlast   = hist['fXaxis']['fNbins'];
      var yfirst  = 1;
      var ylast   = hist['fYaxis']['fNbins'];
      var zfirst  = 1;
      var zlast   = hist['fZaxis']['fNbins'];
      var maximum = Number.NEGATIVE_INFINITY;
      var minimum = Number.POSITIVE_INFINITY;
      var tmp_value;
      for (binz=zfirst;binz<=zlast;binz++) {
         for (biny=yfirst;biny<=ylast;biny++) {
            for (binx=xfirst;binx<=xlast;binx++) {
               //bin = hist.getBin(binx,biny,binz);
               //tmp_value = hist.getBinContent(bin);
               tmp_value = hist.getBinContent(binx, biny);
               if (tmp_value > maximum) maximum = tmp_value;
               if (tmp_value < minimum) minimum = tmp_value;
            }
         }
      }
      hist['fMaximum'] = maximum;
      hist['fMinimum'] = minimum;
      if (what == 'max') return maximum;
      if (what == 'min') return minimum;
   }

   JSROOTPainter.getValueColor = function(hist, zc, options) {
      var wmin = this.getMinMax(hist, 'min'),
          wmax = this.getMinMax(hist, 'max'),
          wlmin = wmin,
          wlmax = wmax,
          ndivz = hist['fContour'].length,
          scale = ndivz / (wlmax - wlmin);
      if (options && options['logz']) {
         if (wmin <= 0 && wmax > 0) wmin = Math.min(1.0, 0.001 * wmax);
         wlmin = Math.log(wmin)/Math.log(10);
         wlmax = Math.log(wmax)/Math.log(10);
      }
      if (default_palette.length == 0) {
         var saturation = 1,
             lightness = 0.5,
             maxHue = 280,
             minHue = 0,
             maxPretty = 50,
             hue;
         for (var i=0 ; i<maxPretty ; i++) {
            hue = (maxHue - (i + 1) * ((maxHue - minHue) / maxPretty))/360.0;
            var rgbval = this.HLStoRGB(hue, lightness, saturation);
            default_palette.push(rgbval);
         }
      }
      if (options && options['logz']) zc = Math.log(zc)/Math.log(10);
      if (zc < wlmin) zc = wlmin;
      var ncolors = default_palette.length;
      var color = Math.round(0.01 + (zc - wlmin) * scale);
      var theColor = Math.round((color + 0.99) * ncolors / ndivz);
      var icol = theColor % ncolors;
      if (icol < 0) icol = 0;
      return default_palette[icol];
   };

   JSROOTPainter.getTimeOffset = function(axis) {

      var timeFormat = axis['fTimeFormat'];
      var i, timeoffset = 0;
      var idF = timeFormat.indexOf('%F');

      if (idF >= 0) {
         timeformat = timeFormat.substr(0, idF);
      } else {
         timeformat = timeFormat;
      }

      if (idF >= 0) {
         var lnF = timeFormat.length;
         var stringtimeoffset = timeFormat.substr(idF+2, lnF);
         for (i=0;i<3;++i) stringtimeoffset = stringtimeoffset.replace('-', '/');
         var stimeoffset = new Date(stringtimeoffset);
         timeoffset = stimeoffset.getTime();
         var ids = stringtimeoffset.indexOf('s');
         if (ids >= 0) {
            var lns = stringtimeoffset.length;
            var sdp = stringtimeoffset.substr(ids+1, lns);
            var dp = parseFloat(sdp);
            timeoffset += dp;
         }
      } else {
         timeoffset = 788918400000; // UTC time at 01/01/95
      }
      return timeoffset;
   };

   JSROOTPainter.formatExp = function(label) {
      var str = label;
      if (parseFloat(str) == 1.0) return '1';
      if (parseFloat(str) == 10.0) return '10';
      var str = str.replace('e+', 'x10@');
      var str = str.replace('e-', 'x10@-');
      var _val = str.substring(0, str.indexOf('@'));
      var _exp = str.substr(str.indexOf('@'));
      _val = _val.replace('@', '');
      _exp = _exp.replace('@', '');
      var u, size = _exp.length;
      for (j=0;j<size;++j) {
         var c = _exp.charAt(j);
         if (c == '+') u = '\u207A';
         else if (c == '-') u = '\u207B'
         else {
            var e = parseInt(c);
            if (e == 1) u = String.fromCharCode(0xB9);
            else if (e > 1 && e < 4) u = String.fromCharCode(0xB0+e);
            else u = String.fromCharCode(0x2070+e);
         }
         _exp = _exp.replace(c, u);
      }
      _val = _val.replace('1x', '');
      return _val+_exp;
   };

   JSROOTPainter.translateExp = function(str) {
      var i, j, lstr = str.match(/\^{[0-9]*}/gi);
      if (lstr != null) {
         var symbol = '';
         for (i=0;i<lstr.length;++i) {
            symbol = lstr[i].replace(' ', '');
            symbol = symbol.replace('^{', ''); // &sup
            symbol = symbol.replace('}', ''); // ;
            var size = symbol.length;
            for (j=0;j<size;++j) {
               var c = symbol.charAt(j);
               var e = parseInt(c);
               if (e == 1) u = String.fromCharCode(0xB9);
               else if (e > 1 && e < 4) u = String.fromCharCode(0xB0+e);
               else u = String.fromCharCode(0x2070+e);
               symbol = symbol.replace(c, u);
            }
            str = str.replace(lstr[i], symbol);
         }
      }
      return str;
   };

   JSROOTPainter.translateLaTeX = function(string) {
      var str = string;
      str = this.translateExp(str);
      while (str.indexOf('^{o}') != -1)
         str = str.replace('^{o}', '\xBA');
      var lstr = str.match(/\#sqrt{(.*?)}/gi);
      if (lstr != null) {
         var symbol;
         for (i=0;i<lstr.length;++i) {
            symbol = lstr[i].replace(' ', '');
            symbol = symbol.replace('#sqrt{', '#sqrt');
            symbol = symbol.replace('}', '');
            str = str.replace(lstr[i], symbol);
         }
      }
      var lstr = str.match(/\_{(.*?)}/gi);
      if (lstr != null) {
         var symbol;
         for (i=0;i<lstr.length;++i) {
            symbol = lstr[i].replace(' ', '');
            symbol = symbol.replace('_{', ''); // &sub
            symbol = symbol.replace('}', ''); // ;
            str = str.replace(lstr[i], symbol);
         }
      }
      var lstr = str.match(/\^{(.*?)}/gi);
      if (lstr != null) {
         var symbol;
         for (i=0;i<lstr.length;++i) {
            symbol = lstr[i].replace(' ', '');
            symbol = symbol.replace('^{', ''); // &sup
            symbol = symbol.replace('}', ''); // ;
            str = str.replace(lstr[i], symbol);
         }
      }
      while (str.indexOf('#/') != -1)
         str = str.replace('#/', symbols_map['#/']);
      for (x in symbols_map) {
         while (str.indexOf(x) != -1)
            str = str.replace(x, symbols_map[x]);
      }
      return str;
   };

   /**
    * Now the real drawing functions (using d3.js)
    */

   JSROOTPainter.add3DInteraction = function(renderer, scene, camera, toplevel) {
      // add 3D mouse interactive functions
      var mouseX, mouseY, mouseDowned = false;
      var mouse = { x: 0, y: 0 }, INTERSECTED;

      var radius = 100;
      var theta = 0;
      var projector = new THREE.Projector();
      function findIntersection() {
         // find intersections
         if ( mouseDowned ) {
            if ( INTERSECTED ) {
               INTERSECTED.material.emissive.setHex( INTERSECTED.currentHex );
               renderer.render( scene, camera );
            }
            INTERSECTED = null;
            tooltip.hide();
            return;
         }
         var vector = new THREE.Vector3( mouse.x, mouse.y, 1 );
         projector.unprojectVector( vector, camera );
         var raycaster = new THREE.Raycaster( camera.position, vector.sub( camera.position ).normalize() );
         var intersects = raycaster.intersectObjects( scene.children, true );
         if ( intersects.length > 0 ) {
            var pick = null;
            for (var i=0;i<intersects.length;++i) {
               if ('emissive' in intersects[ i ].object.material) {
                  pick = intersects[ i ];
                  break;
               }
            }
            if (pick && INTERSECTED != pick.object ) {
               if ( INTERSECTED ) INTERSECTED.material.emissive.setHex( INTERSECTED.currentHex );
               INTERSECTED = pick.object;
               INTERSECTED.currentHex = INTERSECTED.material.emissive.getHex();
               INTERSECTED.material.emissive.setHex( 0x5f5f5f );
               renderer.render( scene, camera );
               tooltip.show(INTERSECTED.name.length > 0 ? INTERSECTED.name : INTERSECTED.parent.name, 200);
            }
         } else {
            if ( INTERSECTED ) {
               INTERSECTED.material.emissive.setHex( INTERSECTED.currentHex );
               renderer.render( scene, camera );
            }
            INTERSECTED = null;
            tooltip.hide();
         }
      };

      $( renderer.domElement ).on('touchstart mousedown',function (e) {
         //var touch = e.changedTouches[0] || {};
         tooltip.hide();
         e.preventDefault();
         var touch = e;
         if ('changedTouches' in e) touch = e.changedTouches[0];
         else if ('touches' in e) touch = e.touches[0];
         else if ('originalEvent' in e) {
            if ('changedTouches' in e.originalEvent) touch = e.originalEvent.changedTouches[0];
            else if ('touches' in e.originalEvent) touch = e.originalEvent.touches[0];
         }
         mouseX = touch.pageX;
         mouseY = touch.pageY;
         mouseDowned = true;
      });
      $( renderer.domElement ).on('touchmove mousemove', function(e) {
         if ( mouseDowned ) {
            var touch = e;
            if ('changedTouches' in e) touch = e.changedTouches[0];
            else if ('touches' in e) touch = e.touches[0];
            else if ('originalEvent' in e) {
               if ('changedTouches' in e.originalEvent) touch = e.originalEvent.changedTouches[0];
               else if ('touches' in e.originalEvent) touch = e.originalEvent.touches[0];
            }
            var moveX = touch.pageX - mouseX;
            var moveY = touch.pageY - mouseY;
            // limited X rotate in -45 to 135 deg
            if ( (moveY > 0 && toplevel.rotation.x < Math.PI*3/4) ||
                 (moveY < 0 &&  toplevel.rotation.x > -Math.PI/4) ) {
               toplevel.rotation.x += moveY*0.02;
            }
            toplevel.rotation.y += moveX*0.02;
            renderer.render( scene, camera );
            mouseX = touch.pageX;
            mouseY = touch.pageY;
         }
         else {
            e.preventDefault();
            var mouse_x = 'offsetX' in e.originalEvent ? e.originalEvent.offsetX : e.originalEvent.layerX;
            var mouse_y = 'offsetY' in e.originalEvent ? e.originalEvent.offsetY : e.originalEvent.layerY;
            mouse.x = ( mouse_x / renderer.domElement.width ) * 2 - 1;
            mouse.y = - ( mouse_y / renderer.domElement.height ) * 2 + 1;
            // enable picking once tootips are available...
            findIntersection();
         }
      });
      $( renderer.domElement ).on('touchend mouseup', function(e) {
         mouseDowned = false;
      });
      $( renderer.domElement ).on('mousewheel', function(e, d) {
         e.preventDefault();
         camera.position.z += d * 20;
         renderer.render( scene, camera );
      });
   };

   JSROOTPainter.addInteraction = function(vis, obj) {
      var width = vis.attr("width"), height = vis.attr("height");
      var e, origin, rect;

      if (!('objects' in vis)) {
         vis['objects'] = new Array();
         doubleTap(vis[0][0]);
      }
      if (vis['objects'].indexOf(obj) != -1) return;
      vis['objects'].push(obj);

      function refresh() {
         if (vis.x_axis && vis.y_axis) {
            vis.select(".xaxis").call(vis.x_axis);
            vis.select(".yaxis").call(vis.y_axis);
         }
         vis.select(".xaxis").selectAll("text")
            .attr("font-size", vis['x_fsize'])
            .attr("font-family", vis['x_font']['name'])
            .attr("font-weight", vis['x_font']['weight'])
            .attr("font-style", vis['x_font']['style']);
         vis.select(".yaxis").selectAll("text")
            .attr("font-size", vis['y_fsize'])
            .attr("font-family", vis['y_font']['name'])
            .attr("font-weight", vis['y_font']['weight'])
            .attr("font-style", vis['y_font']['style']);
         for (var i=0;i<vis['objects'].length;++i) {
            vis['objects'][i].redraw();
         }
      };
      //var zoom = d3.behavior.zoom().x(obj.x).y(obj.y).on("zoom", refresh());
      var zoom = d3.behavior.zoom().x(obj.x).y(obj.y);
      vis.on("touchstart", startRectSel);
      vis.on("mousedown", startRectSel);

      function startRectSel() {
         d3.event.preventDefault();
         vis.select("#zoom_rect").remove();
         e = this;
         var t = d3.event.changedTouches;
         origin = t ? d3.touches(e, t)[0] : d3.mouse(e);
         rect = vis.append("rect").attr("class", "zoom").attr("id", "zoom_rect");
         d3.select("body").classed("noselect", true);
         d3.select("body").style("-webkit-user-select", "none");
         origin[0] = Math.max(0, Math.min(width, origin[0]));
         origin[1] = Math.max(0, Math.min(height, origin[1]));
         vis.on("dblclick", unZoom);
         d3.select(window)
            .on("mousemove.zoomRect", moveRectSel)
            .on("mouseup.zoomRect", endRectSel, true);
         d3.select(window)
            .on("touchmove.zoomRect", moveRectSel)
            .on("touchend.zoomRect", endRectSel, true);
         d3.event.stopPropagation();
      };

      function unZoom() {
         d3.event.preventDefault();
         var xmin = vis['objects'][0]['x_min'],
             xmax = vis['objects'][0]['x_max'],
             ymin = vis['objects'][0]['y_min'],
             ymax = vis['objects'][0]['y_max'];
         for (var i=0;i<vis['objects'].length;++i) {
            zoom.x(vis['objects'][i].x.domain([xmin, xmax]))
                .y(vis['objects'][i].y.domain([ymin, ymax]));
            if ('ys' in vis['objects'][i])
               vis['objects'][i].ys.domain([ymin, ymax])
         }
         refresh();
      };

      function moveRectSel() {
         d3.event.preventDefault();
         var t = d3.event.changedTouches;
         var m = t ? d3.touches(e, t)[0] : d3.mouse(e);
         m[0] = Math.max(0, Math.min(width, m[0]));
         m[1] = Math.max(0, Math.min(height, m[1]));
         rect.attr("x", Math.min(origin[0], m[0]))
             .attr("y", Math.min(origin[1], m[1]))
             .attr("width", Math.abs(m[0] - origin[0]))
             .attr("height", Math.abs(m[1] - origin[1]));
      };

      function endRectSel() {
         d3.event.preventDefault();
         d3.select(window).on("touchmove.zoomRect", null).on("touchend.zoomRect", null);
         d3.select(window).on("mousemove.zoomRect", null).on("mouseup.zoomRect", null);
         d3.select("body").classed("noselect", false);
         var t = d3.event.changedTouches;
         var m = t ? d3.touches(e, t)[0] : d3.mouse(e);
         m[0] = Math.max(0, Math.min(width, m[0]));
         m[1] = Math.max(0, Math.min(height, m[1]));
         if (Math.abs(m[0] - origin[0]) > 10 && Math.abs(m[1] - origin[1]) > 10) {
            var xmin = Math.min(vis['objects'][0].x.invert(origin[0]),
                                vis['objects'][0].x.invert(m[0])),
                xmax = Math.max(vis['objects'][0].x.invert(origin[0]),
                                vis['objects'][0].x.invert(m[0])),
                ymin = Math.min(vis['objects'][0].y.invert(origin[1]),
                                vis['objects'][0].y.invert(m[1])),
                ymax = Math.max(vis['objects'][0].y.invert(origin[1]),
                                vis['objects'][0].y.invert(m[1]));
            for (var i=0;i<vis['objects'].length;++i) {
               zoom.x(vis['objects'][i].x.domain([xmin, xmax]))
                   .y(vis['objects'][i].y.domain([ymin, ymax]));
               if ('ys' in vis['objects'][i])
                  vis['objects'][i].ys.domain([ymin, ymax])
            }
         }
         rect.remove();
         refresh();
         d3.select("body").style("-webkit-user-select", "auto");
      };
   };

   JSROOTPainter.createCanvas = function(element, idx) {
      var w = element.width(), h = w * 0.6666666;
      var render_to = '#histogram' + idx;
      d3.select(render_to).style("background-color", 'white');
      d3.select(render_to).style("width", "100%");

      var svg = d3.select(render_to).append("svg")
                  .attr("width", w)
                  .attr("height", h)
                  .style("background-color", 'white');
      defs = svg.append('svg:defs');
      return svg;
   };

   JSROOTPainter.createFrame = function(vis, pad, histo, frame) {
      var w = vis.attr("width"), h = vis.attr("height");
      var width = w, height = h;
      var lm = w*0.12, rm = w*0.05, tm = h*0.12, bm = h*0.12;
      if (histo && histo['fOption'] && histo['fOption'].toLowerCase() == 'colz')
         rm = w * 0.13;
      var framecolor = 'white', bordermode = 0,
          bordersize = 0, linecolor = 'black',//root_colors[0],
          linestyle = 0, linewidth = 1;
      if (histo && histo['_typename'] == 'JSROOTIO.TF1') {
         linecolor = 'black';
         linewidth = 1;
      }
      if (frame) {
         bordermode = frame['fBorderMode'];
         bordersize = frame['fBorderSize'];
         linecolor = root_colors[frame['fLineColor']];
         linestyle = frame['fLineStyle'];
         linewidth = frame['fLineWidth'];
         if (pad) {
            var xspan = width / Math.abs(pad['fX2'] - pad['fX1']);
            var yspan = height / Math.abs(pad['fY2'] - pad['fY1']);
            px1 = (frame['fX1'] - pad['fX1']) * xspan;
            py1 = (frame['fY1'] - pad['fY1']) * yspan;
            px2 = (frame['fX2'] - pad['fX1']) * xspan;
            py2 = (frame['fY2'] - pad['fY1']) * yspan;
            if (px1 < px2) {pxl = px1; pxt = px2;}
            else           {pxl = px2; pxt = px1;}
            if (py1 < py2) {pyl = py1; pyt = py2;}
            else           {pyl = py2; pyt = py1;}
            lm = pxl;
            bm = pyl;
            w = pxt - pxl;
            h = pyt - pyl;
            tm = height - pyt;
            rm = width - pxt;
         }
         else {
            lm = frame['fX1'] * width;
            tm = frame['fY1'] * height;
            bm = (1.0 - frame['fY2']) * height;
            rm = (1.0 - frame['fX2']) * width;
            w -= (lm + rm);
            h -= (tm + bm);
         }
         framecolor = root_colors[frame['fFillColor']];
         if (frame['fFillStyle'] > 4000 && frame['fFillStyle'] < 4100)
            framecolor = 'none';
      }
      else {
         if (pad) {
            framecolor = root_colors[pad['fFrameFillColor']];
            if (pad['fFrameFillStyle'] > 4000 && pad['fFrameFillStyle'] < 4100)
               framecolor = 'none';
         }
         w -= (lm + rm);
         h -= (tm + bm);
      }
      if (typeof(framecolor) == 'undefined')
         framecolor = 'white';

      var hframe = vis.append("svg:g")
            .attr("x", lm)
            .attr("y", tm)
            .attr("width", w)
            .attr("height", h)
            .attr("transform", "translate(" + lm + "," + tm + ")");

      hframe.append("svg:rect")
            .attr("x", 0)
            .attr("y", 0)
            .attr("width", w)
            .attr("height", h)
            .attr("fill", framecolor)
            .style("stroke", linecolor)
            .style("stroke-width", linewidth);

      var svg_frame = hframe.append("svg")
            .attr("id", "svg_frame_" + (++frame_id))
            .attr("x", 0)
            .attr("y", 0)
            .attr("width", w)
            .attr("height", h)
            .attr("viewBox", "0 0 "+w+" "+h);

      return {
         id: "#svg_frame_" + frame_id,
         frame: hframe,
         xmin: 0,
         xmax: 0,
         ymin: 0,
         ymax: 0
      };
   };

   JSROOTPainter.drawAxes = function(vis, histo, pad, xx, yy) {
      var w = vis.attr("width"), h = vis.attr("height"),
          logx = false, logy = false, logz = false;
      if (pad && typeof(pad) != 'undefined') {
         logx = pad['fLogx'];
         logy = pad['fLogy'];
         logz = pad['fLogz'];
      }
      var noexpx = histo['fXaxis'].TestBit(EAxisBits.kNoExponent);
      var noexpy = histo['fYaxis'].TestBit(EAxisBits.kNoExponent);
      var moreloglabelsx = histo['fXaxis'].TestBit(EAxisBits.kMoreLogLabels);
      var moreloglabelsy = histo['fYaxis'].TestBit(EAxisBits.kMoreLogLabels);

      if (histo['fXaxis']['fXmax'] < 100 && histo['fXaxis']['fXmax']/histo['fXaxis']['fXmin'] < 100) noexpx = true;
      if (histo['fYaxis']['fXmax'] < 100 && histo['fYaxis']['fXmax']/histo['fYaxis']['fXmin'] < 100) noexpy = true;

      var ndivx = histo['fXaxis']['fNdivisions'];
      var n1ax = ndivx%100;
      var n2ax = (ndivx%10000 - n1ax)/100;
      var n3ax = ndivx/10000;

      var ndivy = histo['fYaxis']['fNdivisions'];
      var n1ay = ndivy%100;
      var n2ay = (ndivy%10000 - n1ay)/100;
      var n3ay = ndivy/10000;

      /* X-axis label */
      var label = this.translateLaTeX(histo['fXaxis']['fTitle']);
      var xAxisTitleFontSize = histo['fXaxis']['fTitleSize'] * h;
      var xAxisLabelOffset = 3 + (histo['fXaxis']['fLabelOffset'] * h);
      var xAxisLabelFontSize = histo['fXaxis']['fLabelSize'] * h;
      var xAxisFontDetails = getFontDetails(root_fonts[Math.floor(histo['fXaxis']['fTitleFont']/10)]);

      vis.append("text")
         .attr("class", "X axis label")
         .attr("x", w)
         .attr("y", h)
         .attr("text-anchor", "end")
         .attr("font-family", xAxisFontDetails['name'])
         .attr("font-weight", xAxisFontDetails['weight'])
         .attr("font-style", xAxisFontDetails['style'])
         .attr("font-size", xAxisTitleFontSize)
         .text(label)
         .attr("transform", "translate(0," + (xAxisLabelFontSize + xAxisLabelOffset * histo['fXaxis']['fTitleOffset'] + xAxisTitleFontSize) + ")");

      /* Y-axis label */
      label = this.translateLaTeX(histo['fYaxis']['fTitle']);
      var yAxisTitleFontSize = histo['fYaxis']['fTitleSize'] * h;
      var yAxisLabelOffset = 3 + (histo['fYaxis']['fLabelOffset'] * w);
      var yAxisLabelFontSize = histo['fYaxis']['fLabelSize'] * h;
      var yAxisFontDetails = getFontDetails(root_fonts[Math.floor(histo['fYaxis']['fTitleFont'] /10)]);

      vis.append("text")
         .attr("class", "Y axis label")
         .attr("x", 0)
         .attr("y", -yAxisLabelFontSize - yAxisTitleFontSize - yAxisLabelOffset * histo['fYaxis']['fTitleOffset'])
         .attr("font-family", yAxisFontDetails['name'])
         .attr("font-size", yAxisTitleFontSize)
         .attr("font-weight", yAxisFontDetails['weight'])
         .attr("font-style", yAxisFontDetails['style'])
         .attr("fill", "black")
         .attr("text-anchor", "end")
         .text(label)
         .attr("transform", "rotate(270, 0, 0)");

      var xAxisColor = histo['fXaxis']['fAxisColor'];
      var xDivLength = histo['fXaxis']['fTickLength'] * h;
      var yAxisColor = histo['fYaxis']['fAxisColor'];
      var yDivLength = histo['fYaxis']['fTickLength'] * w;

      /*
       * Define the scales, according to the information from the pad
       */
      var dfx = d3.format(",.f"), dfy = d3.format(",.f");
      if (histo['fXaxis']['fTimeDisplay']) {
         if (n1ax > 8) n1ax = 8;
         var timeoffset = this.getTimeOffset(histo['fXaxis']);
         var range = histo['fXaxis']['fXmax'] - histo['fXaxis']['fXmin'];
         dfx = d3.time.format("%Mm%S");
         if (range>31536000)
            dfx = d3.time.format("%Y");
         else if (range>2419200)
            dfx = d3.time.format("%Y/%m");
         else if (range>86400)
            dfx = d3.time.format("%Y/%m/%d");
         else if (range>3600)
            dfx = d3.time.format("%Hh%Mm%S");
         else if (range>60)
            dfx = d3.time.format("%Hh%M");

         var x_axis = d3.svg.axis()
            .scale(xx)
            .orient("bottom")
            .tickPadding(xAxisLabelOffset)
            .tickSize(-xDivLength, -xDivLength/2, -xDivLength/4)
            .tickFormat(function(d, i) {
               var datime = new Date(timeoffset + (d * 1000));
               return dfx(datime); })
            .ticks(n1ax);
      }
      else if (logx) {
         var x_axis = d3.svg.axis()
            .scale(xx)
            .orient("bottom")
            .tickPadding(xAxisLabelOffset)
            .tickSize(-xDivLength, -xDivLength/2, -xDivLength/4)
            .tickFormat(function(d, i) { var val = parseFloat(d);
               var vlog = Math.abs(JSROOTMath.log10(val));
               if (moreloglabelsx) {
                  if (vlog % 1 < 0.7 || vlog % 1 > 0.9999) {
                     if (noexpx) return val.toFixed();
                     else return JSROOTPainter.formatExp(val.toExponential(0));
                  }
                  else return null;
               }
               else {
                  if (vlog % 1 < 0.0001 || vlog % 1 > 0.9999) {
                     if (noexpx) return val.toFixed();
                     else return JSROOTPainter.formatExp(val.toExponential(0));
                  }
                  else return null;
               }
            });
      }
      else {
         var x_axis = d3.svg.axis()
            .scale(xx)
            .orient("bottom")
            .tickPadding(xAxisLabelOffset)
            .tickSubdivide(n2ax-1)
            .tickSize(-xDivLength, -xDivLength/2, -xDivLength/4)
            .tickFormat(function(d,i) {
               if (histo['fXaxis']['fTimeDisplay']) return dfx;
               return parseFloat(d.toPrecision(12));
            })
            .ticks(n1ax);
      }
      if (histo['fYaxis']['fTimeDisplay']) {
         if (n1ay > 8) n1ay = 8;
         var timeoffset = this.getTimeOffset(histo['fYaxis']);
         var range = histo['fYaxis']['fXmax'] - histo['fYaxis']['fXmin'];
         dfy = d3.time.format("%Mm%S");

         if (range>31536000)
            dfy = d3.time.format("%Y");
         else if (range>2419200)
            dfy = d3.time.format("%Y/%m");
         else if (range>86400)
            dfy = d3.time.format("%Y/%m/%d");
         else if (range>3600)
            dfy = d3.time.format("%Hh%Mm%S");
         else if (range>60)
            dfy = d3.time.format("%Hh%M");

         var y_axis = d3.svg.axis()
            .scale(yy)
            .orient("left")
            .tickPadding(yAxisLabelOffset)
            .tickSize(-yDivLength, -yDivLength/2, -yDivLength/4)
            .tickFormat(function(d, i) {
               var datime = new Date(timeoffset + (d * 1000));
               return dfy(datime); })
            .ticks(n1ay);
      }
      else if (logy) {
         var y_axis = d3.svg.axis()
            .scale(yy)
            .orient("left")
            .tickPadding(yAxisLabelOffset)
            .tickSize(-yDivLength, -yDivLength/2, -yDivLength/4)
            .tickFormat(function(d, i) { var val = parseFloat(d);
               var vlog = Math.abs(JSROOTMath.log10(val));
               if (moreloglabelsy) {
                  if (vlog % 1 < 0.7 || vlog % 1 > 0.9999) {
                     if (noexpy) return val.toFixed();
                     else return JSROOTPainter.formatExp(val.toExponential(0));
                  }
                  else return null;
               }
               else {
                  if (vlog % 1 < 0.0001 || vlog % 1 > 0.9999) {
                     if (noexpy) return val.toFixed();
                     else return JSROOTPainter.formatExp(val.toExponential(0));
                  }
                  else return null;
            }});
      }
      else {
         if (n1ay >= 10) n1ay -= 2;
         var y_axis = d3.svg.axis()
            .scale(yy)
            .orient("left")
            .tickPadding(yAxisLabelOffset)
            .tickSubdivide(n2ay-1)
            .tickSize(-yDivLength, -yDivLength/2, -yDivLength/4)
            .tickFormat(function(d,i) {
               if (histo['fYaxis']['fTimeDisplay']) return dfy;
               return parseFloat(d.toPrecision(12));
            })
            .ticks(n1ay);
      }
      var xax = vis.append("svg:g")
         .attr("class", "xaxis")
         .attr("transform", "translate(0," + h + ")")
         .call(x_axis);

      var yax = vis.append("svg:g")
         .attr("class", "yaxis")
         .call(y_axis);


      var xAxisLabelFontDetails = getFontDetails(root_fonts[Math.floor(histo['fXaxis']['fLabelFont']/10)]);
      var yAxisLabelFontDetails = getFontDetails(root_fonts[Math.floor(histo['fXaxis']['fLabelFont']/10)]);

      xax.selectAll("text")
         .attr("font-family", xAxisLabelFontDetails['name'])
         .attr("font-size", xAxisLabelFontSize)
         .attr("font-weight", xAxisLabelFontDetails['weight'])
         .attr("font-style", xAxisLabelFontDetails['style']);
      yax.selectAll("text")
         .attr("font-family", yAxisLabelFontDetails['name'])
         .attr("font-size", yAxisLabelFontSize)
         .attr("font-weight", yAxisLabelFontDetails['weight'])
         .attr("font-style", yAxisLabelFontDetails['style']);

      vis['x_axis']  = x_axis;
      vis['y_axis']  = y_axis;
      vis['x_fsize'] = xAxisLabelFontSize;
      vis['y_fsize'] = yAxisLabelFontSize;
      vis['x_font'] = xAxisLabelFontDetails;
      vis['y_font'] = yAxisLabelFontDetails;
   };

   JSROOTPainter.drawCanvas = function(canvas, idx) {
      var render_to = '#histogram' + idx,
          w = $(render_to).width(),
          factor = w / Math.abs(canvas['fUtoPixel']),
          h = Math.abs(canvas['fVtoPixel']) * factor,
          fillcolor = root_colors[canvas['fFillColor']];
      if (canvas['fFillStyle'] > 4000 && canvas['fFillStyle'] < 4100)
         fillcolor = 'none';

      d3.select(render_to).style("background-color", fillcolor);
      d3.select(render_to).style("width", "100%");

      var svg = d3.select(render_to).append("svg")
          .attr("width", w)
          .attr("height", h)
          .style("background-color", fillcolor);
      defs = svg.append('svg:defs');

      JSROOTPainter.drawPrimitives(svg, canvas);
      return svg;
   };

   JSROOTPainter.drawErrors = function(svg, bins, histo, pad, x, y) {
      var w = svg.attr("width"), h = svg.attr("height");
      /* Add a panel for each data point */
      var options = JSROOTPainter.decodeOptions(histo['fOption'], histo, pad);
      var info_marker = getRootMarker(root_markers, histo['fMarkerStyle']);
      var shape = info_marker['shape'], filled = info_marker['toFill'],
          toRotate = info_marker['toRotate'], marker_size = histo['fMarkerSize'] * 32;

      if (histo['fMarkerStyle'] == 1) marker_size = 1;

      var marker = d3.svg.symbol()
          .type(d3.svg.symbolTypes[shape])
          .size(marker_size);

      function do_redraw() {

         JSROOTPainter.drawGrid(svg, histo, pad, x, y);

         if (histo['fName'] == '') histo['fName'] = "random_histo_" + random_id++;
         var g_id = format_id(histo['fName']);
         svg.selectAll("#e_"+g_id).remove();
         var g = svg.append("svg:g")
            .attr("id", "e_"+g_id);

         /* Draw x-error indicators */
         g.selectAll("error_x")
            .data(histo.bins)
            .enter()
            .append("svg:line")
            .attr("x1", function(d) { return histo.x(d.x-d.xerr)} )
            .attr("y1", function(d) { return histo.y(d.y)} )
            .attr("x2", function(d) { return histo.x(d.x+d.xerr)} )
            .attr("y2", function(d) { return histo.y(d.y)} )
            .style("stroke", root_colors[histo['fLineColor']])
            .style("stroke-width", histo['fLineWidth']);

         if (options.Error == 11) {
            g.selectAll("e1_x")
               .data(histo.bins)
               .enter()
               .append("svg:line")
               .attr("y1", function(d) { return histo.y(d.y)-3} )
               .attr("x1", function(d) { return histo.x(d.x-d.xerr)})
               .attr("y2", function(d) { return histo.y(d.y)+3})
               .attr("x2", function(d) { return histo.x(d.x-d.xerr)})
               .style("stroke", root_colors[histo['fLineColor']])
               .style("stroke-width", histo['fLineWidth']);
            g.selectAll("e1_x")
               .data(histo.bins)
               .enter()
               .append("svg:line")
               .attr("y1", function(d) { return histo.y(d.y)-3} )
               .attr("x1", function(d) { return histo.x(d.x+d.xerr) })
               .attr("y2", function(d) { return histo.y(d.y)+3})
               .attr("x2", function(d) { return histo.x(d.x+d.xerr) })
               .style("stroke", root_colors[histo['fLineColor']])
               .style("stroke-width", histo['fLineWidth']);
         }
         /* Draw y-error indicators */
         g.selectAll("error_y")
            .data(histo.bins)
            .enter()
            .append("svg:line")
            .attr("x1", function(d) { return histo.x(d.x)})
            .attr("y1", function(d) { return histo.y(d.y-d.yerr) })
            .attr("x2", function(d) { return histo.x(d.x)})
            .attr("y2", function(d) { return histo.y(d.y+d.yerr) })
            .style("stroke", root_colors[histo['fLineColor']])
            .style("stroke-width", histo['fLineWidth']);

         if (options.Error == 11) {
            g.selectAll("e1_y")
               .data(histo.bins)
               .enter()
               .append("svg:line")
               .attr("x1", function(d) { return histo.x(d.x)-3})
               .attr("y1", function(d) { return histo.y(d.y-d.yerr) })
               .attr("x2", function(d) { return histo.x(d.x)+3})
               .attr("y2", function(d) { return histo.y(d.y-d.yerr) })
               .style("stroke", root_colors[histo['fLineColor']])
               .style("stroke-width", histo['fLineWidth']);
            g.selectAll("e1_y")
               .data(histo.bins)
               .enter()
               .append("svg:line")
               .attr("x1", function(d) { return histo.x(d.x)-3})
               .attr("y1", function(d) { return histo.y(d.y+d.yerr) })
               .attr("x2", function(d) { return histo.x(d.x)+3})
               .attr("y2", function(d) { return histo.y(d.y+d.yerr) })
               .style("stroke", root_colors[histo['fLineColor']])
               .style("stroke-width", histo['fLineWidth']);
         }
         var points = g.selectAll("markers")
            .data(histo.bins)
            .enter()
            .append("svg:path")
            .attr("class", "marker")
            .attr("transform", function(d) {
               return "translate(" + histo.x(d.x) + "," + histo.y(d.y) + ")"
            })
            .style("fill", root_colors[histo['fMarkerColor']])
            .style("stroke", root_colors[histo['fMarkerColor']])
            .attr("d", marker)
            .append("svg:title")
            .text(function(d) {
               return "x = " + d.x.toPrecision(4) + " \ny = " + d.y.toPrecision(4) +
                      " \nerror x = " + d.xerr.toPrecision(4) +
                      " \nerror y = " + d.yerr.toPrecision(4);
            });
         g.selectAll("line")
           .append("svg:title")
           .text(function(d) {
              return "x = " + d.x.toPrecision(4) + " \ny = " + d.y.toPrecision(4) +
                     " \nerror x = " + d.xerr.toPrecision(4) +
                     " \nerror y = " + d.yerr.toPrecision(4);
           });
      };
      histo['redraw'] = do_redraw;
      do_redraw();
   };

   JSROOTPainter.drawFunction = function(vis, pad, func, hframe) {
      var i, logx = false, logy = false, logz = false,
          gridx = false, gridy = false, draw_all = true;
      if (pad && typeof(pad) != 'undefined') {
         logx = pad['fLogx'];
         logy = pad['fLogy'];
         logz = pad['fLogz'];
         gridx = pad['fGridx'];
         gridy = pad['fGridy'];
      }
      var fillcolor = root_colors[func['fFillColor']];
      var linecolor = root_colors[func['fLineColor']];
      if (func['fFillColor'] == 0) {
         fillcolor = '#4572A7';
      }
      if (func['fLineColor'] == 0) {
         linecolor = '#4572A7';
      }
      var interpolate_method = 'basis';
      var h, hmin = 1.0e32, hmax = -1.0e32;
      if (func['fNsave'] > 0) {
         // in the case where the points have been saved, useful for example
         // if we don't have the user's function
         var nb_points = func['fNpx'];
         for (var i=0;i<nb_points;++i) {
            h = func['fSave'][i];
            if (h > hmax) hmax = h;
            if (h < hmin) hmin = h;
         }
         if (hmax > 0.0) hmax *= 1.05;
         if (hmin < 0.0) hmin *= 1.05;
         func['fYmin'] = hmin;
         func['fYmax'] = hmax;
         func['x_min'] = func['fSave'][nb_points+1];
         func['x_max'] = func['fSave'][nb_points+2];
         func['y_min'] = func['fYmin'];
         func['y_max'] = func['fYmax'];
         binwidth = ((func['x_max'] - func['x_min']) / nb_points);
         var bins = d3.range(nb_points).map(function(p) {
            return {
               x: func['x_min'] + (p * binwidth),
               y: func['fSave'][p]
            };
         });
         func['bins'] = bins;
         interpolate_method = 'monotone';
         //interpolate_method = 'cardinal-open';
      }
      else {
         // we don't have the points, so let's try to interpret the function
         // use fNpfits instead of fNpx if possible (to use more points)
         if (func['fNpfits'] <= 103) func['fNpfits'] = 333;
         var nb_points = Math.max(func['fNpx'], func['fNpfits']);
         var binwidth = ((func['fXmax'] - func['fXmin']) / nb_points);
         for (var i=0;i<nb_points;++i) {
            h = func.evalPar(func['fXmin'] + (i * binwidth));
            if (isNaN(h)) {
               func['fXmin'] = 1.2e-308;
               h = func.evalPar(func['fXmin'] + (i * binwidth));
            }
            if (h > hmax) hmax = h;
            if (h < hmin) hmin = h;
         }
         if (hmax > 0.0) hmax *= 1.05;
         if (hmin < 0.0) hmin *= 1.05;
         func['fYmin'] = hmin;
         func['fYmax'] = hmax;
         func['x_min'] = func['fXmin'];
         func['x_max'] = func['fXmax'];
         func['y_min'] = func['fYmin'];
         func['y_max'] = func['fYmax'];
         var bins = d3.range(nb_points).map(function(p) {
            return {
               x: func['fXmin'] + (p * binwidth),
               y: func.evalPar(func['fXmin'] + (p * binwidth))
            };
         });
         func['bins'] = bins;
         interpolate_method = 'cardinal-open';
      }
      var ret = hframe != null ? hframe : this.createFrame(vis, pad, func, null);
      var frame = ret['frame'];
      var svg_frame = d3.select(ret['id']);
      var w = frame.attr("width"), h = frame.attr("height");
      if (hframe == null || (hframe['xmin'] < 1e-300 && hframe['xmax'] < 1e-300 &&
          hframe['ymin'] < 1e-300 && hframe['ymax'] < 1e-300)) {
         if (logx)
            var x = d3.scale.log().domain([func['fXmin'], func['fXmax']]).range([0, w]);
         else
            var x = d3.scale.linear().domain([func['fXmin'], func['fXmax']]).range([0, w]);
         if (logy)
            var y = d3.scale.log().domain([hmin, hmax]).range([h, 0]);
         else
            var y = d3.scale.linear().domain([hmin, hmax]).range([h, 0]);
      }
      else {
         draw_all = false;
         if (logx)
            var x = d3.scale.log().domain([hframe['xmin'], hframe['xmax']]).range([0, w]);
         else
            var x = d3.scale.linear().domain([hframe['xmin'], hframe['xmax']]).range([0, w]);
         if (logy)
            var y = d3.scale.log().domain([hframe['ymin'], hframe['ymax']]).range([h, 0]);
         else
            var y = d3.scale.linear().domain([hframe['ymin'], hframe['ymax']]).range([h, 0]);
      }
      func['x'] = x;
      func['y'] = y;

      function do_redraw() {

         if (func['fName'] == '') func['fName'] = "random_function_" + random_id++;
         var g_id = format_id(func['fName']);
         svg_frame.selectAll("#"+g_id).remove();

         var g = svg_frame.append("svg:g")
            .attr("id", g_id);

         var line = d3.svg.line()
            .x(function(d) { return func.x(d.x);})
            .y(function(d) { return func.y(d.y);})
            .interpolate(interpolate_method);

         g.append("svg:path")
            .attr("class", "line")
            .attr("d", line(func.bins))
            .style("stroke", linecolor)
            .style("stroke-width", func['fLineWidth'])
            .style("stroke-dasharray", root_line_styles[func['fLineStyle']])
            .style("fill", "none");

         // add tooltips
         g.selectAll("line")
            .data(bins)
            .enter()
            .append("svg:circle")
            .attr("cx", function(d) { return x(d.x); })
            .attr("cy", function(d) { return y(d.y); })
            .attr("r", 3)
            .attr("opacity", 0)
            .append("svg:title").text(function(d) {
               return "x = " + d.x.toPrecision(4) + " \ny = " + d.y.toPrecision(4);
            });
      };
      func['redraw'] = do_redraw;
      do_redraw();

      if (draw_all) {

         var mul = (w > h) ? h : w;
         var label_font_size = Math.round(0.035 * mul);

         /* X-axis */
         var x_axis = d3.svg.axis()
            .scale(x)
            .orient("bottom")
            .tickPadding(5)
            .tickSubdivide(4)
            .tickSize(-0.03 * h, -0.03 * h / 2, null)
            .tickFormat(function(d,i) { return parseFloat(d.toPrecision(12)); })
            .ticks(10);

         /* Y-axis minor ticks */
         var y_axis = d3.svg.axis()
            .scale(y)
            .orient("left")
            .tickPadding(5)
            .tickSubdivide(4)
            .tickSize(-0.03 * w, -0.03 * w / 2, null)
            .tickFormat(function(d,i) { return parseFloat(d.toPrecision(12)); })
            .ticks(8);

         var xax = frame.append("svg:g")
            .attr("class", "xaxis")
            .attr("transform", "translate(0," + h + ")")
            .call(x_axis);

         var yax = frame.append("svg:g")
            .attr("class", "yaxis")
            .call(y_axis);

         var font_size = Math.round(0.050 * h);

         if (!pad || typeof(pad) == 'undefined') {
            vis.append("text")
               .attr("class", "title")
               .attr("text-anchor", "middle")
               .attr("x", vis.attr("width")/2)
               .attr("y", 1 + font_size)
               .attr("font-family", "Arial")
               .attr("font-size", font_size)
               .text(func['fTitle']);
         }
         xax.selectAll("text").attr("font-size", label_font_size);
         yax.selectAll("text").attr("font-size", label_font_size);

         frame['x_axis']  = x_axis;
         frame['y_axis']  = y_axis;
         frame['x_fsize'] = label_font_size;
         frame['y_fsize'] = label_font_size;
         frame['x_font']  = {'weight' : "",'style' : "", 'name' : "arial" };
         frame['y_font']  = {'weight' : "",'style' : "", 'name' : "arial" };
      }
      this.addInteraction(frame, func);
      func_list.push(func);
   };

   JSROOTPainter.drawFunctions = function(vis, histo, pad, frame) {
      /* draw statistics box & other TPaveTexts */
      var draw_stats = true;
      if ('fFunctions' in histo) {
         for (i=0; i<histo['fFunctions'].length; ++i) {
            if (histo['fFunctions'][i]['_typename'] == 'JSROOTIO.TPaveStats')
               draw_stats = false;
            if (histo['fFunctions'][i]['_typename'] == 'JSROOTIO.TPaveText' ||
                histo['fFunctions'][i]['_typename'] == 'JSROOTIO.TPaveStats') {
               if (histo['fFunctions'][i]['fX1NDC'] < 1.0 && histo['fFunctions'][i]['fY1NDC'] < 1.0 &&
                   histo['fFunctions'][i]['fX1NDC'] > 0.0 && histo['fFunctions'][i]['fY1NDC'] > 0.0) {
                  this.drawPaveText(vis, histo['fFunctions'][i]);
               }
            }
            if (histo['fFunctions'][i]['_typename'] == 'JSROOTIO.TF1') {
               if (!pad && !histo['fFunctions'][i].TestBit(kNotDraw)) {
                  //if (histo['fFunctions'][i].TestBit(EStatusBits.kObjInCanvas)) {
                     if (typeof(histo['fFunctions'][i]['isDrawn']) == 'undefined' ||
                         histo['fFunctions'][i]['isDrawn'] == false)
                        this.drawFunction(vis, pad, histo['fFunctions'][i], frame);
                     histo['fFunctions'][i]['isDrawn'] = true;
                  //}
               }
               else if (pad && histo['fFunctions'][i].TestBit(EStatusBits.kObjInCanvas)) {
                  if (typeof(histo['fFunctions'][i]['isDrawn']) == 'undefined' ||
                      histo['fFunctions'][i]['isDrawn'] == false)
                     this.drawFunction(vis, pad, histo['fFunctions'][i], frame);
                  histo['fFunctions'][i]['isDrawn'] = true;
               }
            }
         }
      }
      return draw_stats;
   };

   JSROOTPainter.drawGraph = function(vis, pad, graph, hframe) {
      var logx = false, logy = false, logz = false,
          gridx = false, gridy = false, draw_all = true;
      var optionLine, optionAxis, optionCurve, optionStar, optionMark,
          optionBar, optionR, optionOne, optionE, optionFill, optionZ,
          optionCurveFill;
      var draw_errors = true;
      if ('fOption' in graph) {
         var opt = graph['fOption'].toUpperCase();
         opt.replace('SAME', '');
      }
      else var opt = 'LP';

      if (opt.indexOf('L') != -1) optionLine = 1;  else optionLine = 0;
      if (opt.indexOf('A') != -1) optionAxis = 1;  else optionAxis = 0;
      if (opt.indexOf('C') != -1) optionCurve= 1;  else optionCurve= 0;
      if (opt.indexOf('*') != -1) optionStar = 1;  else optionStar = 0;
      if (opt.indexOf('P') != -1) optionMark = 1;  else optionMark = 0;
      if (opt.indexOf('B') != -1) optionBar  = 1;  else optionBar  = 0;
      if (opt.indexOf('R') != -1) optionR    = 1;  else optionR    = 0;
      if (opt.indexOf('1') != -1) optionOne  = 1;  else optionOne  = 0;
      if (opt.indexOf('F') != -1) optionFill = 1;  else optionFill = 0;
      if (opt.indexOf('2') != -1 || opt.indexOf('3') != -1 ||
          opt.indexOf('4') != -1 || opt.indexOf('5') != -1) optionE = 1;
      else optionE = 0;
      optionZ = 0;

      // if no drawing option is selected and if chopt<>' ' nothing is done.
      if (optionLine + optionFill + optionCurve + optionStar + optionMark + optionBar + optionE == 0) {
         if (chopt.length == 0) optionLine = 1;
         else return;
      }
      if (optionStar) graph['fMarkerStyle'] = 3;
      optionCurveFill = 0;
      if (optionCurve && optionFill) {
         optionCurveFill = 1;
         optionFill      = 0;
      }
      if (pad && typeof(pad) != 'undefined') {
         logx = pad['fLogx'];
         logy = pad['fLogy'];
         logz = pad['fLogz'];
         gridx = pad['fGridx'];
         gridy = pad['fGridy'];
      }
      var xaxis_type = logx ? 'logarithmic' : 'linear';
      if (graph['_typename'] == 'JSROOTIO.TGraph') {
         // check for axis scale format, and convert if required
         if (graph['fHistogram']['fXaxis']['fTimeDisplay']) {
            xaxis_type = 'datetime';
         }
         var yaxis_type = logy ? 'logarithmic' : 'linear';
         if (graph['fHistogram']['fYaxis']['fTimeDisplay']) {
            yaxis_type = 'datetime';
         }
      }
      else if (graph['_typename'] == 'JSROOTIO.TGraphErrors') {
         maxEX = d3.max(graph['fEX']);
         maxEY = d3.max(graph['fEY']);
         if (maxEX < 1.0e-300 && maxEY < 1.0e-300)
            draw_errors = false;
      }
      var seriesType = 'scatter';
      if (optionBar == 1)
         seriesType = 'bar';
      var showMarker = false;
      if (optionMark == 1 || optionStar == 1)
         showMarker = true;
      if (optionLine == 1 || optionCurve == 1 || optionFill == 1)
         seriesType = 'line';

      if (optionBar == 1) {
         var binwidth = (graph['fHistogram']['fXaxis']['fXmax'] - graph['fHistogram']['fXaxis']['fXmin']) /
                         graph['fNpoints'];
      }
      var bins = d3.range(graph['fNpoints']).map(function(p) {
         if (optionBar == 1) {
            return {
               x: graph['fX'][p] - (binwidth / 2),
               y: graph['fY'][p], // graph['fHistogram']['fXaxis']['fXmin'],
               bw: binwidth,
               bh: graph['fY'][p]
            }
         }
         else if (graph['_typename'] == 'JSROOTIO.TGraphErrors') {
            return {
               x: graph['fX'][p],
               y: graph['fY'][p],
               exlow:  graph['fEX'][p],
               exhigh: graph['fEX'][p],
               eylow:  graph['fEY'][p],
               eyhigh: graph['fEY'][p]
            };
         }
         else if (graph['_typename'] == 'JSROOTIO.TGraphAsymmErrors' ||
             graph['_typename'].match(/\bRooHist/)) {
            return {
               x: graph['fX'][p],
               y: graph['fY'][p],
               exlow:  graph['fEXlow'][p],
               exhigh: graph['fEXhigh'][p],
               eylow:  graph['fEYlow'][p],
               eyhigh: graph['fEYhigh'][p]
            };
         }
         else {
            return {
               x: graph['fX'][p],
               y: graph['fY'][p]
            };
         }
      });
      var ret = hframe != null ? hframe : this.createFrame(vis, pad, graph['fHistogram'], null);
      var frame = ret['frame'];
      var svg_frame = d3.select(ret['id']);
      var w = frame.attr("width"), h = frame.attr("height");
      if (hframe == null || (hframe['xmin'] < 1e-300 && hframe['xmax'] < 1e-300 &&
          hframe['ymin'] < 1e-300 && hframe['ymax'] < 1e-300)) {
         if (logx)
            var x = d3.scale.log().domain([graph['fHistogram']['fXaxis']['fXmin'],
                                 graph['fHistogram']['fXaxis']['fXmax']]).range([0, w]);
         else
            var x = d3.scale.linear().domain([graph['fHistogram']['fXaxis']['fXmin'],
                                    graph['fHistogram']['fXaxis']['fXmax']]).range([0, w]);
         if (logy)
            var y = d3.scale.log().domain([graph['fHistogram']['fYaxis']['fXmin'],
                                 graph['fHistogram']['fYaxis']['fXmax']]).range([h, 0]);
         else
            var y = d3.scale.linear().domain([graph['fHistogram']['fYaxis']['fXmin'],
                                    graph['fHistogram']['fYaxis']['fXmax']]).range([h, 0]);
         if (logy)
            var ys = d3.scale.log().domain([graph['fHistogram']['fYaxis']['fXmin'],
                                 graph['fHistogram']['fYaxis']['fXmax']]).range([0, h]);
         else
            var ys = d3.scale.linear().domain([graph['fHistogram']['fYaxis']['fXmin'],
                                    graph['fHistogram']['fYaxis']['fXmax']]).range([0, h]);
         graph['x_min'] = graph['fHistogram']['fXaxis']['fXmin'];
         graph['x_max'] = graph['fHistogram']['fXaxis']['fXmax'];
         graph['y_min'] = graph['fHistogram']['fYaxis']['fXmin'];
         graph['y_max'] = graph['fHistogram']['fYaxis']['fXmax'];
      }
      else {
         draw_all = false;
         if (logx)
            var x = d3.scale.log().domain([hframe['xmin'], hframe['xmax']]).range([0, w]);
         else
            var x = d3.scale.linear().domain([hframe['xmin'], hframe['xmax']]).range([0, w]);
         if (logy)
            var y = d3.scale.log().domain([hframe['ymin'], hframe['ymax']]).range([h, 0]);
         else
            var y = d3.scale.linear().domain([hframe['ymin'], hframe['ymax']]).range([h, 0]);

         if (logy)
            var ys = d3.scale.log().domain([hframe['ymin'], hframe['ymax']]).range([0, h]);
         else
            var ys = d3.scale.linear().domain([hframe['ymin'], hframe['ymax']]).range([0, h]);

         graph['x_min'] = hframe['xmin'];
         graph['x_max'] = hframe['xmax'];
         graph['y_min'] = hframe['ymin'];
         graph['y_max'] = hframe['ymax'];
      }
      graph['x'] = x;
      graph['y'] = y;
      graph['ys'] = ys;
      graph['bins'] = bins;

      // exclusion graphs
      var lw = graph['fLineWidth'];
      var ec, ff = 1, exclusionGraph = false;
      if (graph['fLineWidth'] > 99) {
         exclusionGraph = true;
         var normx, normy;
         var n = graph['fNpoints'];
         var glw = graph['fLineWidth'],
             xo = new Array(n+2), yo = new Array(n+2),
             xt = new Array(n+2), yt = new Array(n+2),
             xf = new Array(2*n+2), yf = new Array(2*n+2);
         // negative value means another side of the line...
         if (glw > 32767) {
            glw = 65536 - glw;
         }
         lw = glw % 100; // line width
         if (lw > 0) optionLine = 1;
         ec = root_colors[graph['fFillColor']];
         ec = ec.replace('rgb', 'rgba');
         ec = ec.replace(')', ', 0.20)');

         var a, i, j, nf, wk = (glw/100)*0.005;
         if (graph['fLineWidth'] > 32767)
            wk *= -1;

         var ratio = w / h;

         var xmin = graph['x_min'], xmax = graph['x_max'],
             ymin = graph['y_min'], ymax = graph['y_max'];
         for (i=0; i<n; i++) {
            xo[i] = (graph['fX'][i] - xmin) / (xmax - xmin);
            yo[i] = (graph['fY'][i] - ymin) / (ymax - ymin);
            if (w > h) yo[i] = yo[i] / ratio;
            else if (h > w) xo[i] = xo[i] / ratio;
         }
         // The first part of the filled area is made of the graph points.
         // Make sure that two adjacent points are different.
         xf[0] = xo[0];
         yf[0] = yo[0];
         nf = 0;
         for (i=1; i<n; i++) {
            if (xo[i] == xo[i-1] && yo[i] == yo[i-1]) continue;
            nf++;
            xf[nf] = xo[i];
            if (xf[i] == xf[i-1]) xf[i] += 0.000001; // add an epsilon to avoid exact vertical lines.
            yf[nf] = yo[i];
         }
         // For each graph points a shifted points is computed to build up
         // the second part of the filled area. First and last points are
         // treated as special cases, outside of the loop.
         if (xf[1] == xf[0]) {
            a = Math.PI / 2.0;
         } else {
            a = Math.atan((yf[1] - yf[0]) / (xf[1] - xf[0]));
         }
         if (xf[0] <= xf[1]) {
            xt[0] = xf[0] - wk * Math.sin(a);
            yt[0] = yf[0] + wk * Math.cos(a);
         } else {
            xt[0] = xf[0] + wk * Math.sin(a);
            yt[0] = yf[0] - wk * Math.cos(a);
         }
         if (xf[nf] == xf[nf-1]) {
            a = Math.PI / 2.0;
         } else {
            a = Math.atan((yf[nf] - yf[nf-1]) / (xf[nf] - xf[nf-1]));
         }
         if (xf[nf] >= xf[nf-1]) {
            xt[nf] = xf[nf] - wk * Math.sin(a);
            yt[nf] = yf[nf] + wk * Math.cos(a);
         } else {
            xt[nf] = xf[nf] + wk * Math.sin(a);
            yt[nf] = yf[nf] - wk * Math.cos(a);
         }

         var a1, a2, a3, xi0, yi0, xi1, yi1, xi2, yi2;
         for (i=1; i<nf; i++) {
            xi0 = xf[i];
            yi0 = yf[i];
            xi1 = xf[i+1];
            yi1 = yf[i+1];
            xi2 = xf[i-1];
            yi2 = yf[i-1];
            if (xi1 == xi0) {
               a1 = Math.PI / 2.0;
            } else {
               a1  = Math.atan((yi1 - yi0) / (xi1 - xi0));
            }
            if (xi1 < xi0) a1 = a1 + Math.PI;
            if (xi2 == xi0) {
               a2 = Math.PI / 2.0;
            } else {
               a2  = Math.atan((yi0 - yi2) / (xi0 - xi2));
            }
            if (xi0 < xi2) a2 = a2 + Math.PI;
            x1 = xi0 - wk * Math.sin(a1);
            y1 = yi0 + wk * Math.cos(a1);
            x2 = xi0 - wk * Math.sin(a2);
            y2 = yi0 + wk * Math.cos(a2);
            xm = (x1 + x2) * 0.5;
            ym = (y1 + y2) * 0.5;
            if (xm == xi0) {
               a3 = Math.PI / 2.0;
            } else {
               a3 = Math.atan((ym - yi0) / (xm - xi0));
            }
            x3 = xi0 - wk * Math.sin(a3 + (Math.PI / 2.0));
            y3 = yi0 + wk * Math.cos(a3 + (Math.PI / 2.0));
            // Rotate (x3,y3) by PI around (xi0,yi0) if it is not on the (xm,ym) side.
            if ((xm - xi0) * (x3 - xi0) < 0 && (ym - yi0) * (y3 - yi0) < 0) {
               x3 = 2 * xi0 - x3;
               y3 = 2 * yi0 - y3;
            }
            if ((xm == x1) && (ym == y1)) {
               x3 = xm;
               y3 = ym;
            }
            xt[i] = x3;
            yt[i] = y3;
         }
         // Close the polygon if the first and last points are the same
         if (xf[nf] == xf[0] && yf[nf] == yf[0]) {
            xm = (xt[nf] + xt[0]) * 0.5;
            ym = (yt[nf] + yt[0]) * 0.5;
            if (xm == xf[0]) {
               a3 = Math.PI / 2.0;
            } else {
               a3 = Math.atan((ym - yf[0]) / (xm - xf[0]));
            }
            x3 = xf[0] + wk * Math.sin(a3 + (Math.PI / 2.0));
            y3 = yf[0] - wk * Math.cos(a3 + (Math.PI / 2.0));
            if ((xm - xf[0]) * (x3 - xf[0]) < 0 && (ym - yf[0]) * (y3 - yf[0]) < 0) {
               x3 = 2 * xf[0] - x3;
               y3 = 2 * yf[0] - y3;
            }
            xt[nf] = x3;
            xt[0]  = x3;
            yt[nf] = y3;
            yt[0]  = y3;
         }
         // Find the crossing segments and remove the useless ones
         var xc, yc, c1, b1, c2, b2;
         var cross = false;
         var nf2 = nf;
         for (i=nf2; i>0; i--) {
            for (j=i-1; j>0; j--) {
               if (xt[i-1] == xt[i] || xt[j-1] == xt[j]) continue;
               c1  = (yt[i-1] - yt[i]) / (xt[i-1] - xt[i]);
               b1  = yt[i] - c1 * xt[i];
               c2  = (yt[j-1] - yt[j]) / (xt[j-1] - xt[j]);
               b2  = yt[j] - c2 * xt[j];
               if (c1 != c2) {
                  xc = (b2 - b1) / (c1 - c2);
                  yc = c1 * xc + b1;
                  if (xc > Math.min(xt[i], xt[i-1]) && xc < Math.max(xt[i], xt[i-1]) &&
                      xc > Math.min(xt[j], xt[j-1]) && xc < Math.max(xt[j], xt[j-1]) &&
                      yc > Math.min(yt[i], yt[i-1]) && yc < Math.max(yt[i], yt[i-1]) &&
                      yc > Math.min(yt[j], yt[j-1]) && yc < Math.max(yt[j], yt[j-1])) {
                     nf++; xf[nf] = xt[i]; yf[nf] = yt[i];
                     nf++; xf[nf] = xc   ; yf[nf] = yc;
                     i = j;
                     cross = true;
                     break;
                  } else {
                     continue;
                  }
               } else {
                  continue;
               }
            }
            if (!cross) {
               nf++;
               xf[nf] = xt[i];
               yf[nf] = yt[i];
            }
            cross = false;
         }
         nf++; xf[nf] = xt[0]; yf[nf] = yt[0]; nf++;
         for (i=0; i<nf; i++) {
            if (w > h) {
               xf[i] = xmin + (xf[i] * (xmax - xmin));
               yf[i] = ymin + (yf[i] * (ymax - ymin)) * ratio;
            }
            else if (h > w) {
               xf[i] = xmin + (xf[i] * (xmax - xmin)) * ratio;
               yf[i] = ymin + (yf[i] * (ymax - ymin));
            }
            else {
               xf[i] = xmin + (xf[i] * (xmax - xmin));
               yf[i] = ymin + (yf[i] * (ymax - ymin));
            }
            if (logx && xf[i] <= 0.0) xf[i] = xmin;
            if (logy && yf[i] <= 0.0) yf[i] = ymin;
         }
         var excl = d3.range(nf).map(function(p) {
            return {
               x: xf[p],
               y: yf[p]
            };
         });
         /* some clean-up */
         xo.splice(0, xo.length); yo.splice(0, yo.length);
         xo = null; yo = null;
         xt.splice(0, xt.length); yt.splice(0, yt.length);
         xt = null; yt = null;
         xf.splice(0, xf.length); yf.splice(0, yf.length);
         xf = null; yf = null;
      }

      function do_redraw() {

         if (draw_all)
            JSROOTPainter.drawGrid(frame, graph['fHistogram'], pad, x, y);

         if (graph['fName'] == '') graph['fName'] = "random_graph_" + random_id++;
         var g_id = format_id(graph['fName']);
         svg_frame.selectAll("#"+g_id).remove();
         var g = svg_frame.append("svg:g")
            .attr("id", g_id);

         if (seriesType == 'line') {
            /* contour lines only */
            var line = d3.svg.line()
               .x(function(d) { return graph.x(d.x);})
               .y(function(d) { return graph.y(d.y);});
         }
         if (seriesType == 'bar') {
            var fillcolor = root_colors[graph['fFillColor']];
            if (typeof(fillcolor) == 'undefined') fillcolor = "rgb(204, 204, 204)";
            /* filled bar graph */
            var xdom = graph.x.domain();
            var xfactor = xdom[1]-xdom[0];
            g.selectAll("bar_graph")
               .data(graph.bins)
               .enter()
               .append("svg:rect")
               .attr("x", function(d) { return graph.x(d.x)} )
               .attr("y", function(d) { return graph.y(d.y)} )
               .attr("width", function(d) { return (w / (xdom[1]-xdom[0]))-1} )
               .attr("height", function(d) { return graph.ys(d.bh)} )
               .style("fill", fillcolor)
               .append("svg:title").text(function(d) {
                  return "x = " + d.x.toPrecision(4) + " \nentries = " + d.y.toPrecision(4);
               });
         }
         if (exclusionGraph) {
            /* first draw exclusion area, and then the line */
            showMarker = false;
            if (graph['fFillStyle'] > 3000 && graph['fFillStyle'] <= 3025) {
               createFillPatterns(vis, graph['fFillStyle'], graph['fFillColor']);
               g.append("svg:path")
                  .attr("class", "line")
                  .attr("d", line(excl))
                  .style("stroke", "none")
                  .style("stroke-width", ff)
                  .style("fill", "url(#pat" + graph['fFillStyle'] + "_" + graph['fFillColor'] + ")")
                  .style("antialias", "false");
            }
            else {
               g.append("svg:path")
                  .attr("class", "line")
                  .attr("d", line(excl))
                  .style("stroke", "none")
                  .style("stroke-width", ff)
                  .style("fill", ec);
            }
         }
         if (seriesType == 'line') {
            g.append("svg:path")
               .attr("class", "line")
               .attr("d", line(bins))
               .style("stroke", (optionLine == 1) ? root_colors[graph['fLineColor']] : "none")
               .style("stroke-width", lw)
               .style("stroke-dasharray", root_line_styles[graph['fLineStyle']])
               .style("fill", (optionFill == 1) ? root_colors[graph['fFillColor']] : "none");

            // add tooltips
            g.selectAll("line")
               .data(bins)
               .enter()
               .append("svg:circle")
               .attr("cx", function(d) { return x(d.x); })
               .attr("cy", function(d) { return y(d.y); })
               .attr("r", 3)
               .attr("opacity", 0)
               .append("svg:title").text(function(d) {
                  return "x = " + d.x.toPrecision(4) + " \ny = " + d.y.toPrecision(4);
               });
         }
         if ((graph['_typename'] == 'JSROOTIO.TGraphErrors' ||
              graph['_typename'] == 'JSROOTIO.TGraphAsymmErrors' ||
              graph['_typename'].match(/\bRooHist/)) && draw_errors && !optionBar) {
            /* Add x-error indicators */
            g.selectAll("error_x")
               .data(graph.bins)
               .enter()
               .append("svg:line")
               .attr("x1", function(d) { return graph.x(d.x-d.exlow)} )
               .attr("y1", function(d) { return graph.y(d.y)} )
               .attr("x2", function(d) { return graph.x(d.x+d.exhigh)} )
               .attr("y2", function(d) { return graph.y(d.y)} )
               .style("stroke", root_colors[graph['fLineColor']])
               .style("stroke-width", graph['fLineWidth']);

            g.selectAll("e1_x")
               .data(graph.bins)
               .enter()
               .append("svg:line")
               .attr("y1", function(d) { return graph.y(d.y)-3} )
               .attr("x1", function(d) { return graph.x(d.x-d.exlow)})
               .attr("y2", function(d) { return graph.y(d.y)+3})
               .attr("x2", function(d) { return graph.x(d.x-d.exlow)})
               .style("stroke", root_colors[graph['fLineColor']])
               .style("stroke-width", graph['fLineWidth']);
            g.selectAll("e1_x")
               .data(graph.bins)
               .enter()
               .append("svg:line")
               .attr("y1", function(d) { return graph.y(d.y)-3} )
               .attr("x1", function(d) { return graph.x(d.x+d.exhigh) })
               .attr("y2", function(d) { return graph.y(d.y)+3})
               .attr("x2", function(d) { return graph.x(d.x+d.exhigh) })
               .style("stroke", root_colors[graph['fLineColor']])
               .style("stroke-width", graph['fLineWidth']);

            /* Add y-error indicators */
            g.selectAll("error_y")
               .data(graph.bins)
               .enter()
               .append("svg:line")
               .attr("x1", function(d) { return graph.x(d.x)})
               .attr("y1", function(d) { return graph.y(d.y-d.eylow) })
               .attr("x2", function(d) { return graph.x(d.x)})
               .attr("y2", function(d) { return graph.y(d.y+d.eyhigh) })
               .style("stroke", root_colors[graph['fLineColor']])
               .style("stroke-width", graph['fLineWidth']);

            g.selectAll("e1_y")
               .data(graph.bins)
               .enter()
               .append("svg:line")
               .attr("x1", function(d) { return graph.x(d.x)-3})
               .attr("y1", function(d) { return graph.y(d.y-d.eylow) })
               .attr("x2", function(d) { return graph.x(d.x)+3})
               .attr("y2", function(d) { return graph.y(d.y-d.eylow) })
               .style("stroke", root_colors[graph['fLineColor']])
               .style("stroke-width", graph['fLineWidth']);
            g.selectAll("e1_y")
               .data(graph.bins)
               .enter()
               .append("svg:line")
               .attr("x1", function(d) { return graph.x(d.x)-3})
               .attr("y1", function(d) { return graph.y(d.y+d.eyhigh) })
               .attr("x2", function(d) { return graph.x(d.x)+3})
               .attr("y2", function(d) { return graph.y(d.y+d.eyhigh) })
               .style("stroke", root_colors[graph['fLineColor']])
               .style("stroke-width", graph['fLineWidth']);

            g.selectAll("line")
              .append("svg:title")
              .text(function(d) {
                  return "x = " + d.x.toPrecision(4) + " \ny = " + d.y.toPrecision(4) +
                         "\nerror x = -" + d.exlow.toPrecision(4) + "/+" + d.exhigh.toPrecision(4) +
                         "\nerror y = -" + d.eylow.toPrecision(4) + "/+" + d.eyhigh.toPrecision(4);
               });
         }
         else draw_errors = false;
         if (showMarker) {
            /* Add markers */
            var filled = false;
            if ((graph['fMarkerStyle'] == 8) ||
                (graph['fMarkerStyle'] > 19 && graph['fMarkerStyle'] < 24) ||
                (graph['fMarkerStyle'] == 29))
               filled = true;

            var info_marker = getRootMarker(root_markers, graph['fMarkerStyle']);

            var shape = info_marker['shape'];
            var filled = info_marker['toFill'];
            var toRotate = info_marker['toRotate'];
            var markerSize = graph['fMarkerSize'];
            var markerScale = (shape == 0) ? 32 : 64;
            if (graph['fMarkerStyle'] == 1) markerScale = 1;

            switch (shape) {
               case 6:
                  var marker = "M " + (-4 * markerSize) + " " + (-1 * markerSize)
                              + " L " + 4 * markerSize + " " + (-1 * markerSize)
                              + " L " + (-2.4 * markerSize) + " " + 4 * markerSize
                              + " L 0 " + (-4 * markerSize) + " L " + 2.8 * markerSize
                              + " " + 4 * markerSize + " z";
                  break;
               case 7:
                  var marker = "M " + (- 4 * markerSize) + " " + (-4 * markerSize)
                              + " L " + 4 * markerSize + " " + 4 * markerSize + " M 0 "
                              + (-4 * markerSize) + " 0 " + 4 * markerSize + " M "
                              + 4 * markerSize + " " + (-4 * markerSize) + " L "
                              + (-4 * markerSize) + " " + 4 * markerSize + " M "
                              + (-4 * markerSize) + " 0 L " + 4 * markerSize + " 0";
                  break;
               default:
                  var marker = d3.svg.symbol()
                              .type(d3.svg.symbolTypes[shape])
                              .size(markerSize * markerScale);
                  break;
            }
            g.selectAll("markers")
               .data(graph.bins)
               .enter()
               .append("svg:path")
               .attr("class", "marker")
               .attr("transform", function(d) {return "translate(" + graph.x(d.x) + "," + graph.y(d.y) + ")"})
               .style("fill", filled ? root_colors[graph['fMarkerColor']] : "none")
               .style("stroke", root_colors[graph['fMarkerColor']])
               .attr("d", marker)
               .append("svg:title")
               .text(function(d) {
                  if (draw_errors)
                     return "x = " + d.x.toPrecision(4) + " \ny = " + d.y.toPrecision(4) +
                            "\nerror x = -" + d.exlow.toPrecision(4) + "/+" + d.exhigh.toPrecision(4) +
                            "\nerror y = -" + d.eylow.toPrecision(4) + "/+" + d.eyhigh.toPrecision(4);
                  else
                     return "x = " + d.x.toPrecision(4) + " \ny = " + d.y.toPrecision(4);
               });

         }
      };
      graph['redraw'] = do_redraw;
      do_redraw();

      if (draw_all) {
         this.drawAxes(frame, graph['fHistogram'], pad, x, y);
         this.drawTitle(vis, graph['fHistogram'], pad);
      }
      this.addInteraction(frame, graph);
      if ('fHistogram' in graph && graph['fHistogram'])
         this.drawFunctions(vis, graph['fHistogram'], pad, ret);
   };

   JSROOTPainter.drawGrid = function(vis, histo, pad, x, y) {
      var gridx = false, gridy = false;
      if (pad && typeof(pad) != 'undefined') {
         gridx = pad['fGridx'];
         gridy = pad['fGridy'];
      }
      var ndivx = histo['fXaxis']['fNdivisions'];
      var n1ax = ndivx%100;
      var n2ax = (ndivx%10000 - n1ax)/100;
      var n3ax = ndivx/10000;
      var nn3x = Math.max(n3ax,1);
      var nn2x = Math.max(n2ax,1)*nn3x;
      var nn1x = Math.max(n1ax,1)*nn2x;

      var ndivy = histo['fYaxis']['fNdivisions'];
      var n1ay = ndivy%100;
      var n2ay = (ndivy%10000 - n1ay)/100;
      var n3ay = ndivy/10000;
      var nn3y = Math.max(n3ay,1);
      var nn2y = Math.max(n2ay,1)*nn3y;
      var nn1y = Math.max(n1ay,1)*nn2y;

      vis.selectAll(".gridLine").remove();

      /* add a grid on x axis, if the option is set */
      if (gridx) {
         vis.selectAll("gridLine")
            .data(x.ticks(n1ax))
            .enter()
            .append("svg:line")
            .attr("class", "gridLine")
            .attr("x1", x)
            .attr("y1", vis.attr("height"))
            .attr("x2", x)
            .attr("y2", 0)
            .style("stroke", "black")
            .style("stroke-width", histo['fLineWidth'])
            .style("stroke-dasharray", root_line_styles[11]);
      }

      /* add a grid on y axis, if the option is set */
      if (gridy) {
         vis.selectAll("gridLine")
            .data(y.ticks(n1ay))
            .enter()
            .append("svg:line")
            .attr("class", "gridLine")
            .attr("x1", 0)
            .attr("y1", y)
            .attr("x2", vis.attr("width"))
            .attr("y2", y)
            .style("stroke", "black")
            .style("stroke-width", histo['fLineWidth'])
            .style("stroke-dasharray", root_line_styles[11]);
      }
   };

   JSROOTPainter.drawHistogram1D = function(vis, pad, histo, hframe) {
      var i, gridx = false, gridy = false;
      var options = JSROOTPainter.decodeOptions(histo['fOption'], histo, pad);
      var draw_all = false;
      if (hframe == null || (hframe['xmin'] < 1e-300 && hframe['xmax'] < 1e-300 &&
          hframe['ymin'] < 1e-300 && hframe['ymax'] < 1e-300)) {
         draw_all = true;
      }
      if (pad && typeof(pad) != 'undefined') {
         gridx = pad['fGridx'];
         gridy = pad['fGridy'];
      }
      var fillcolor = root_colors[histo['fFillColor']];
      var linecolor = root_colors[histo['fLineColor']];
      if (histo['fFillColor'] == 0) {
         fillcolor = '#4572A7';
      }
      if (histo['fLineColor'] == 0) {
         linecolor = '#4572A7';
      }
      var hmin = 1.0e32, hmax = -1.0e32;
      for (i=0;i<histo['fXaxis']['fNbins'];++i) {
         if (histo['fArray'][i+1] < hmin) hmin = histo['fArray'][i+1];
         if (histo['fArray'][i+1] > hmax) hmax = histo['fArray'][i+1];
      }
      var mul = (hmin < 0) ? 1.05 : 1.0;
      if (Math.abs(hmin) < 1e-300 && Math.abs(hmax) < 1e-300) {
         var ymin = histo['fYaxis']['fXmin'], ymax = histo['fYaxis']['fXmax'];
         if (histo['fMinimum'] != -1111) ymin = histo['fMinimum'];
         if (histo['fMaximum'] != -1111) ymax = histo['fMaximum'];

         // special case used for drawing multiple graphs in the same frame
         var ret = hframe != null ? hframe : this.createFrame(vis, pad, histo, null);
         var frame = ret['frame'];
         var svg_frame = d3.select(ret['id']);
         var w = frame.attr("width"), h = frame.attr("height");
         if (options.Logx)
            var x = d3.scale.log().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([0, w]);
         else
            var x = d3.scale.linear().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([0, w]);
         if (options.Logy)
            var y = d3.scale.log().domain([ymin, ymax]).range([h, 0]);
         else
            var y = d3.scale.linear().domain([ymin, ymax]).range([h, 0]);

         // avoid this!
         histo['fYaxis']['fXmin'] = ymin;
         histo['fYaxis']['fXmax'] = ymax;

         histo['x'] = x;
         histo['y'] = y;
         histo['x_min'] = histo['fXaxis']['fXmin'];
         histo['x_max'] = histo['fXaxis']['fXmax'];
         histo['y_min'] = ymin;
         histo['y_max'] = ymax;
         histo['redraw'] = function() {
            JSROOTPainter.drawGrid(frame, histo, pad, x, y);
         };
         this.drawAxes(frame, histo, pad, x, y);
         this.drawTitle(vis, histo, pad);
         this.addInteraction(frame, histo);
         this.drawFunctions(vis, histo, pad, ret);
         histo.redraw();
         return {
            frame: frame,
            xmin: histo['fXaxis']['fXmin'],
            xmax: histo['fXaxis']['fXmax'],
            ymin: ymin,
            ymax: ymax
         };
      }
      if (histo['fMinimum'] != -1111) hmin = histo['fMinimum'];
      if (histo['fMaximum'] != -1111) hmax = histo['fMaximum'];
      histo['fYaxis']['fXmin'] = hmin * mul;
      histo['fYaxis']['fXmax'] = hmax * 1.05;
      var binwidth = ((histo['fXaxis']['fXmax'] - histo['fXaxis']['fXmin']) / histo['fXaxis']['fNbins']);
      var bins = d3.range(histo['fXaxis']['fNbins']).map(function(p) {
         var offset = (options.Error > 0) ? (p * binwidth) - (binwidth / 2.0) : (p * binwidth);
         return {
            x:  histo['fXaxis']['fXmin'] + offset,
            y:  histo['fArray'][p],
            xerr: binwidth / 2.0,
            yerr: histo.getBinError(p)
         };
      });
      var ret = hframe != null ? hframe : this.createFrame(vis, pad, histo, null);
      var frame = ret['frame'];
      var svg_frame = d3.select(ret['id']);
      var w = frame.attr("width"), h = frame.attr("height");
      if (options.Logx)
         var x = d3.scale.log().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([0, w]);
      else
         var x = d3.scale.linear().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([0, w]);
      if (options.Logy)
         var y = d3.scale.log().domain([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]).range([h, 0]);
      else
         var y = d3.scale.linear().domain([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]).range([h, 0]);

      if (options.Same) {
         x.domain([ret['xmin'],ret['xmax']]);
         y.domain([ret['ymin'],ret['ymax']]);
      }
      else {
         ret['xmin'] = histo['fXaxis']['fXmin'],
         ret['xmax'] = histo['fXaxis']['fXmax'],
         ret['ymin'] = histo['fYaxis']['fXmin'];
         ret['ymax'] = histo['fYaxis']['fXmax'];
      }
      if (histo['fXaxis'].TestBit(EAxisBits.kAxisRange)) {
         ret['xmin'] = histo.getBinLowEdge(histo['fXaxis']['fFirst']);
         ret['xmax'] = histo.getBinUpEdge(histo['fXaxis']['fLast']);
         x.domain([ret['xmin'],ret['xmax']]);
         y.domain([ret['ymin'],ret['ymax']]);
      }
      histo['x_min'] = histo['fXaxis']['fXmin'];
      histo['x_max'] = histo['fXaxis']['fXmax'];
      histo['y_min'] = histo['fYaxis']['fXmin'];
      histo['y_max'] = histo['fYaxis']['fXmax'];

      histo['x'] = x;
      histo['y'] = y;
      histo['bins'] = bins;

      if (options.Error > 0) {
         this.drawErrors(svg_frame, bins, histo, pad, x, y);
      }
      else if (options.Bar == 0 && options.Hist == 0) {
         histo['redraw'] = function() {
            JSROOTPainter.drawGrid(frame, histo, pad, x, y);
         };
      }
      else {
         function do_redraw() {

            if (draw_all)
               JSROOTPainter.drawGrid(frame, histo, pad, x, y);

            if (histo['fName'] == '') histo['fName'] = "random_histo_" + random_id++;
            var g_id = format_id(histo['fName']);
            svg_frame.selectAll("#"+g_id).remove();
            var g = svg_frame.append("svg:g")
               .attr("id", g_id);

            if ((histo['fFillStyle'] < 4000 || histo['fFillStyle'] > 4100) && histo['fFillColor'] != 0) {

               /* histogram filling */
               var area = d3.svg.area()
                  .x(function(d) { return histo.x(d.x);})
                  .y0(function(d) { return histo.y(0);})
                  .y1(function(d) { return histo.y(d.y);})
                  .interpolate("step-before")

               if (histo['fFillStyle'] > 3000 && histo['fFillStyle'] <= 3025) {
                  createFillPatterns(vis, histo['fFillStyle'], histo['fFillColor']);
                  g.append("svg:path")
                     .attr("class", "area")
                     .attr("d", area(bins))
                     .style("stroke", linecolor)
                     .style("stroke-width", histo['fLineWidth'])
                     .style("fill", "url(#pat" + histo['fFillStyle'] + "_" + histo['fFillColor'] + ")")
                     .style("antialias", "false");
               }
               else {
                  g.append("svg:path")
                     .attr("class", "area")
                     .attr("d", area(bins))
                     .style("stroke", linecolor)
                     .style("stroke-width", histo['fLineWidth'])
                     .style("fill", fillcolor)
                     .style("antialias", "false");
               }
            }
            else {
               /* histogram contour lines only */
               var line = d3.svg.line()
                  .x(function(d) { return histo.x(d.x);})
                  .y(function(d) { return histo.y(d.y);})
                  .interpolate("step-before");

               g.append("svg:path")
                  .attr("class", "line")
                  .attr("d", line(bins))
                  .style("stroke", linecolor)
                  .style("stroke-width", histo['fLineWidth'])
                  .style("fill", "none")
                  .style("stroke-dasharray", histo['fLineStyle'] > 1 ? root_line_styles[histo['fLineStyle']] : null)
                  .style("antialias", "false");
            }
            // add tooltips
            var selwidth = x(2*binwidth)-x(binwidth);
            g.selectAll("selections")
               .data(bins)
               .enter()
               .append("svg:line")
               .attr("x1", function(d) { return x(d.x-d.xerr) } )
               .attr("y1", function(d) { return y(d.y) } )
               .attr("x2", function(d) { return x(d.x-d.xerr) } )
               .attr("y2", function(d) { return y(0) } )
               .attr("opacity", 0)
               .style("stroke", "#4572A7")
               .style("stroke-width", selwidth)
               .on('mouseover', function() { d3.select(this).transition().duration(100).style("opacity", 0.3) } )
               .on('mouseout', function() { d3.select(this).transition().duration(100).style("opacity", 0) } )
               .append("svg:title").text(function(d) { return "x = [" + (d.x-(2*d.xerr)).toPrecision(4) +
                       ", " + d.x.toPrecision(4) + "] \nentries = " + d.y;
               });
         };
         histo['redraw'] = do_redraw;
         do_redraw();
      }
      if (draw_all)
         this.drawAxes(frame, histo, pad, x, y);
      this.drawTitle(vis, histo, pad);
      this.addInteraction(frame, histo);
      var draw_stats = this.drawFunctions(vis, histo, pad, ret);
      if (draw_stats && !pad || typeof(pad) == 'undefined')
         this.drawStat(vis, histo);
      return null;
   };

   JSROOTPainter.drawHistogram2D = function(vis, pad, histo, hframe) {
      var i, gridx = false, gridy = false;
      var options = JSROOTPainter.decodeOptions(histo['fOption'], histo, pad);
      if (pad && typeof(pad) != 'undefined') {
         gridx = pad['fGridx'];
         gridy = pad['fGridy'];
      }
      var fillcolor = root_colors[histo['fFillColor']];
      var linecolor = root_colors[histo['fLineColor']];
      if (histo['fFillColor'] == 0) {
         fillcolor = '#4572A7';
      }
      if (histo['fLineColor'] == 0) {
         linecolor = '#4572A7';
      }
      var nbinsx = histo['fXaxis']['fNbins'];
      var nbinsy = histo['fYaxis']['fNbins'];
      var scalex = (histo['fXaxis']['fXmax'] - histo['fXaxis']['fXmin']) /
                    histo['fXaxis']['fNbins'];
      var scaley = (histo['fYaxis']['fXmax'] - histo['fYaxis']['fXmin']) /
                    histo['fYaxis']['fNbins'];
      var maxbin = -1e32, minbin = 1e32;
      for (i=0; i<nbinsx; ++i) {
         for (j=0; j<nbinsy; ++j) {
            var bin_content = histo.getBinContent(i, j);
            if (bin_content < minbin) minbin = bin_content;
            if (bin_content > maxbin) maxbin = bin_content;
         }
      }
      var bins = new Array();
      for (i=0; i<nbinsx; ++i) {
         for (var j=0; j<nbinsy; ++j) {
            var bin_content = histo.getBinContent(i, j);
            if (bin_content > minbin) {
               var point = {
                  x:histo['fXaxis']['fXmin'] + (i*scalex),
                  y:histo['fYaxis']['fXmin'] + (j*scaley),
                  z:bin_content
               };
               bins.push(point);
            }
         }
      }
      var ret = hframe != null ? hframe : this.createFrame(vis, pad, histo, null);
      var frame = ret['frame'];
      var svg_frame = d3.select(ret['id']);
      var w = frame.attr("width"), h = frame.attr("height");
      if (options.Logx)
         var x = d3.scale.log().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([0, w]);
      else
         var x = d3.scale.linear().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([0, w]);
      if (options.Logy)
         var y = d3.scale.log().domain([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]).range([h, 0]);
      else
         var y = d3.scale.linear().domain([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]).range([h, 0]);

      var c = d3.scale.linear().domain([minbin, maxbin]).range(['red', 'blue']);

      histo['x_min'] = histo['fXaxis']['fXmin'];
      histo['x_max'] = histo['fXaxis']['fXmax'];
      histo['y_min'] = histo['fYaxis']['fXmin'];
      histo['y_max'] = histo['fYaxis']['fXmax'];

      histo['x'] = x;
      histo['y'] = y;
      histo['bins'] = bins;

      function do_redraw() {

         JSROOTPainter.drawGrid(frame, histo, pad, x, y);

         if (histo['fName'] == '') histo['fName'] = "random_histo_" + random_id++;
         var g_id = format_id(histo['fName']);
         svg_frame.selectAll("#"+g_id).remove();
         var g = svg_frame.append("svg:g")
            .attr("id", g_id);

         var xdom = histo.x.domain();
         var ydom = histo.y.domain();
         var xfactor = Math.abs(histo['fXaxis']['fXmax']-histo['fXaxis']['fXmin']) / Math.abs(xdom[1]-xdom[0]);
         var yfactor = Math.abs(histo['fYaxis']['fXmax']-histo['fYaxis']['fXmin']) / Math.abs(ydom[1]-ydom[0]);

         if (options.Scat > 0 && histo['fMarkerStyle'] > 1) {
            /* Add markers */
            var filled = false;
            if ((histo['fMarkerStyle'] == 8) ||
                (histo['fMarkerStyle'] > 19 && histo['fMarkerStyle'] < 24) ||
                (histo['fMarkerStyle'] == 29))
               filled = true;

            var info_marker = getRootMarker(root_markers, histo['fMarkerStyle']);

            var shape = info_marker['shape'];
            var filled = info_marker['toFill'];
            var toRotate = info_marker['toRotate'];
            var markerSize = histo['fMarkerSize'];
            var markerScale = (shape == 0) ? 32 : 64;
            if (histo['fMarkerStyle'] == 1) markerScale = 1;

            switch (shape) {
               case 6:
                  var marker = "M " + (-4 * markerSize) + " " + (-1 * markerSize)
                              + " L " + 4 * markerSize + " " + (-1 * markerSize)
                              + " L " + (-2.4 * markerSize) + " " + 4 * markerSize
                              + " L 0 " + (-4 * markerSize) + " L " + 2.8 * markerSize
                              + " " + 4 * markerSize + " z";
                  break;
               case 7:
                  var marker = "M " + (- 4 * markerSize) + " " + (-4 * markerSize)
                              + " L " + 4 * markerSize + " " + 4 * markerSize + " M 0 "
                              + (-4 * markerSize) + " 0 " + 4 * markerSize + " M "
                              + 4 * markerSize + " " + (-4 * markerSize) + " L "
                              + (-4 * markerSize) + " " + 4 * markerSize + " M "
                              + (-4 * markerSize) + " 0 L " + 4 * markerSize + " 0";
                  break;
               default:
                  var marker = d3.svg.symbol()
                              .type(d3.svg.symbolTypes[shape])
                              .size(markerSize * markerScale);
                  break;
            }
            g.selectAll("markers")
               .data(histo.bins)
               .enter()
               .append("svg:path")
               .attr("class", "marker")
               .attr("transform", function(d) {
                  return "translate(" + histo.x(d.x) + "," + histo.y(d.y) + ")"
               })
               .style("fill", root_colors[histo['fMarkerColor']])
               .style("stroke", root_colors[histo['fMarkerColor']])
               .attr("d", marker)
               .append("svg:title")
               .text(function(d) { return "x = " + d.x.toPrecision(4) + " \ny = " +
                                   d.y.toPrecision(4) + " \nentries = " + d.z; });
         }
         else {
            g.selectAll("bins")
               .data(histo.bins)
               .enter()
               .append("svg:rect")
               .attr("class", "bins")
               .attr("x", function(d) {
                  return histo.x(d.x + (scalex / 2)) - (0.5 * d.z * ((w / histo['fXaxis']['fNbins']) / maxbin) * xfactor);
               })
               .attr("y", function(d) {
                  return histo.y(d.y + (scaley / 2)) - (0.5 * d.z * ((h / histo['fYaxis']['fNbins']) / maxbin) * yfactor);
               })
               .attr("width", function(d) {
                  if (options.Color > 0)
                     return (w / histo['fXaxis']['fNbins']) * xfactor;
                  else
                     return d.z * ((w / histo['fXaxis']['fNbins']) / maxbin) * xfactor;
               })
               .attr("height", function(d) {
                  if (options.Color > 0)
                     return (h / histo['fYaxis']['fNbins']) * yfactor;
                  else
                     return d.z * ((h / histo['fYaxis']['fNbins']) / maxbin) * yfactor;
               })
               .style("stroke", function(d) {
                  if (options.Color > 0)
                     return JSROOTPainter.getValueColor(histo, d.z, pad);
                  else
                     return "black";
               })
               .style("fill", function(d) {
                  if (options.Color > 0)
                     return JSROOTPainter.getValueColor(histo, d.z, pad);
                  else
                     return "none";
               });
            g.selectAll("selections")
               .data(bins)
               .enter()
               .append("svg:rect")
               .attr("x", function(d) {
                  return histo.x(d.x + (scalex / 2)) - (0.5 * d.z * ((w / histo['fXaxis']['fNbins']) / maxbin) * xfactor);
               })
               .attr("y", function(d) {
                  return histo.y(d.y + (scaley / 2)) - (0.5 * d.z * ((h / histo['fYaxis']['fNbins']) / maxbin) * yfactor);
               })
               .attr("width", function(d) {
                  if (options.Color > 0)
                     return (w / histo['fXaxis']['fNbins']) * xfactor;
                  else
                     return d.z * ((w / histo['fXaxis']['fNbins']) / maxbin) * xfactor;
               })
               .attr("height", function(d) {
                  if (options.Color > 0)
                     return (h / histo['fYaxis']['fNbins']) * yfactor;
                  else
                     return d.z * ((h / histo['fYaxis']['fNbins']) / maxbin) * yfactor;
               })
               .attr("opacity", 0)
               .style("stroke", "#4572A7")
               .style("fill", "#4572A7")
               .on('mouseover', function() { d3.select(this).transition().duration(100).style("opacity", 0.3) } )
               .on('mouseout', function() { d3.select(this).transition().duration(100).style("opacity", 0) } )
               .append("svg:title")
               .text(function(d) {
                  return "x = [" + d.x.toPrecision(4) + ", " + (d.x + scalex).toPrecision(4) +
                  "] \ny = [" + d.y.toPrecision(4) + ", " + (d.y + scaley).toPrecision(4) +
                  "] \nentries = " + d.z;
                });
         }
      };
      histo['redraw'] = do_redraw;
      do_redraw();

      if (options.Color > 0 && options.Zscale > 0) {
         // just to initialize the default palette
         this.getValueColor(histo, 0, pad);
         for (i=0; i<histo['fFunctions'].length; ++i) {
            if (histo['fFunctions'][i]['_typename'] == 'JSROOTIO.TPaletteAxis')
               this.drawPaletteAxis(vis, histo['fFunctions'][i], minbin, maxbin);
         }
      }

      this.drawAxes(frame, histo, pad, x, y);
      this.drawTitle(vis, histo, pad);
      this.addInteraction(frame, histo);
      var draw_stats = this.drawFunctions(vis, histo, pad, ret);
      if (draw_stats && !pad || typeof(pad) == 'undefined')
         this.drawStat(vis, histo);
   };

   JSROOTPainter.drawHistogram2D3D = function(vis, pad, histo, hframe) {
      var i, j, k, logx = false, logy = false, logz = false,
          gridx = false, gridy = false, girdz = false;
      var opt = histo['fOption'].toLowerCase();
      if (pad && typeof(pad) != 'undefined') {
         logx = pad['fLogx'];
         logy = pad['fLogy'];
         logz = pad['fLogz'];
         gridx = pad['fGridx'];
         gridy = pad['fGridy'];
         gridz = pad['fGridz'];
      }
      var fillcolor = root_colors[histo['fFillColor']];
      var linecolor = root_colors[histo['fLineColor']];
      if (histo['fFillColor'] == 0) {
         fillcolor = '#4572A7';
      }
      if (histo['fLineColor'] == 0) {
         linecolor = '#4572A7';
      }
      var nbinsx = histo['fXaxis']['fNbins'];
      var nbinsy = histo['fYaxis']['fNbins'];
      var scalex = (histo['fXaxis']['fXmax'] - histo['fXaxis']['fXmin']) /
                    histo['fXaxis']['fNbins'];
      var scaley = (histo['fYaxis']['fXmax'] - histo['fYaxis']['fXmin']) /
                    histo['fYaxis']['fNbins'];
      var maxbin = -1e32, minbin = 1e32;
      for (i=0; i<nbinsx; ++i) {
         for (j=0; j<nbinsy; ++j) {
            var bin_content = histo.getBinContent(i, j);
            if (bin_content < minbin) minbin = bin_content;
            if (bin_content > maxbin) maxbin = bin_content;
         }
      }
      maxbin *= 1.05;
      var bins = new Array();
      for (i=0; i<nbinsx; ++i) {
         for (j=0; j<nbinsy; ++j) {
            var bin_content = histo.getBinContent(i, j);
            if (bin_content > minbin) {
               var point = {
                  x:histo['fXaxis']['fXmin'] + (i*scalex),
                  y:histo['fYaxis']['fXmin'] + (j*scaley),
                  z:bin_content
               };
               bins.push(point);
            }
         }
      }
      var w = vis.attr("width"), h = vis.attr("height"), size = 100;
      if (logx) {
         var tx = d3.scale.log().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([-size, size]);
         var utx = d3.scale.log().domain([-size, size]).range([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]);
      } else {
         var tx = d3.scale.linear().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([-size, size]);
         var utx = d3.scale.linear().domain([-size, size]).range([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]);
      }
      if (logy) {
         var ty = d3.scale.log().domain([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]).range([-size, size]);
         var uty = d3.scale.log().domain([size, -size]).range([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]);
      } else {
         var ty = d3.scale.linear().domain([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]).range([-size, size]);
         var uty = d3.scale.linear().domain([size, -size]).range([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]);
      }
      if (logz) {
         var tz = d3.scale.log().domain([minbin, Math.ceil( maxbin/10 )*10]).range([0, size*2]);
         var utz = d3.scale.log().domain([0, size*2]).range([minbin, Math.ceil( maxbin/10 )*10]);
      } else {
         var tz = d3.scale.linear().domain([minbin, Math.ceil( maxbin/10 )*10]).range([0, size*2]);
         var utz = d3.scale.linear().domain([0, size*2]).range([minbin, Math.ceil( maxbin/10 )*10]);
      }

      // three.js 3D drawing
      var scene = new THREE.Scene();

      var toplevel = new THREE.Object3D();
      toplevel.rotation.x = 30 * Math.PI / 180;
      toplevel.rotation.y = 30 * Math.PI / 180;
      scene.add( toplevel );

      var wireMaterial = new THREE.MeshBasicMaterial( {
         color: 0x000000,
         wireframe: true,
         wireframeLinewidth: 0.5,
         side: THREE.DoubleSide } );

      // create a new mesh with cube geometry
      var cube = new THREE.Mesh( new THREE.CubeGeometry( size*2, size*2, size*2 ), wireMaterial);
      cube.position.y = size;

      // add the cube to the scene
      toplevel.add( cube );

      var textMaterial = new THREE.MeshBasicMaterial( { color: 0x000000 } );

      // add the calibration vectors and texts
      var geometry = new THREE.Geometry();
      var imax, istep, len = 3, plen, sin45 = Math.sin(45);
      var text3d, text;
      var xmajors = tx.ticks(8);
      var xminors = tx.ticks(50);
      for ( i=-size, j=0, k=0; i<size; ++i ) {
         var is_major = ( utx( i ) <= xmajors[j] && utx( i+1 ) > xmajors[j] ) ? true : false;
         var is_minor = ( utx( i ) <= xminors[k] && utx( i+1 ) > xminors[k] ) ? true : false;
         plen = ( is_major ? len + 2 : len) * sin45;
         if ( is_major ) {
            text3d = new THREE.TextGeometry( xmajors[j], {
               size: 7,
               height: 0,
               curveSegments: 10
            });
            ++j;

            text3d.computeBoundingBox();
            var centerOffset = 0.5 * ( text3d.boundingBox.max.x - text3d.boundingBox.min.x );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( i-centerOffset, -13, size+plen );
            toplevel.add( text );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( i+centerOffset, -13, -size-plen );
            text.rotation.y = Math.PI;
            toplevel.add( text );
         }
         if ( is_major || is_minor ) {
            ++k;
            geometry.vertices.push( new THREE.Vector3( i, 0, size ) );
            geometry.vertices.push( new THREE.Vector3( i, -plen, size+plen ) );
            geometry.vertices.push( new THREE.Vector3( i, 0, -size ) );
            geometry.vertices.push( new THREE.Vector3( i, -plen, -size-plen ) );
         }
      }
      var ymajors = ty.ticks(8);
      var yminors = ty.ticks(50);
      for ( i=size, j=0, k=0; i>-size; --i ) {
         var is_major = ( uty( i ) <= ymajors[j] && uty( i-1 ) > ymajors[j] ) ? true : false;
         var is_minor = ( uty( i ) <= yminors[k] && uty( i-1 ) > yminors[k] ) ? true : false;
         plen = ( is_major ? len + 2 : len) * sin45;
         if ( is_major ) {
            text3d = new THREE.TextGeometry( ymajors[j], {
               size: 7,
               height: 0,
               curveSegments: 10
            });
            ++j;

            text3d.computeBoundingBox();
            var centerOffset = 0.5 * ( text3d.boundingBox.max.x - text3d.boundingBox.min.x );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( size+plen, -13, i+centerOffset );
            text.rotation.y = Math.PI/2;
            toplevel.add( text );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( -size-plen, -13, i-centerOffset );
            text.rotation.y = -Math.PI/2;
            toplevel.add( text );
         }
         if ( is_major || is_minor ) {
            ++k;
            geometry.vertices.push( new THREE.Vector3( size, 0, i ) );
            geometry.vertices.push( new THREE.Vector3( size+plen, -plen, i ) );
            geometry.vertices.push( new THREE.Vector3( -size, 0, i ) );
            geometry.vertices.push( new THREE.Vector3( -size-plen, -plen, i ) );
         }
      }
      var zmajors = tz.ticks(8);
      var zminors = tz.ticks(50);
      for ( i=0, j=0, k=0; i<(size*2); ++i ) {
         var is_major = ( utz( i ) <= zmajors[j] && utz( i+1 ) > zmajors[j] ) ? true : false;
         var is_minor = ( utz( i ) <= zminors[k] && utz( i+1 ) > zminors[k] ) ? true : false;
         plen = ( is_major ? len + 2 : len) * sin45;
         if ( is_major ) {
            text3d = new THREE.TextGeometry( zmajors[j], {
               size: 7,
               height: 0,
               curveSegments: 10
            });
            ++j;

            text3d.computeBoundingBox();
            var offset = 0.8 * ( text3d.boundingBox.max.x - text3d.boundingBox.min.x );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( size+offset+5, i-2.5, size+offset+5 );
            text.rotation.y = Math.PI*3/4;
            toplevel.add( text );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( size+offset+5, i-2.5, -size-offset-5 );
            text.rotation.y = -Math.PI*3/4;
            toplevel.add( text );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( -size-offset-5, i-2.5, size+offset+5 );
            text.rotation.y = Math.PI/4;
            toplevel.add( text );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( -size-offset-5, i-2.5, -size-offset-5 );
            text.rotation.y = -Math.PI/4;
            toplevel.add( text );
         }
         if ( is_major || is_minor ) {
            ++k;
            geometry.vertices.push( new THREE.Vector3( size, i, size ) );
            geometry.vertices.push( new THREE.Vector3( size+plen, i, size+plen ) );
            geometry.vertices.push( new THREE.Vector3( size, i, -size ) );
            geometry.vertices.push( new THREE.Vector3( size+plen, i, -size-plen ) );
            geometry.vertices.push( new THREE.Vector3( -size, i, size ) );
            geometry.vertices.push( new THREE.Vector3( -size-plen, i, size+plen ) );
            geometry.vertices.push( new THREE.Vector3( -size, i, -size ) );
            geometry.vertices.push( new THREE.Vector3( -size-plen, i, -size-plen ) );
         }
      }

      // add the calibration lines
      var lineMaterial = new THREE.LineBasicMaterial( { color: 0x000000 } );
      var line = new THREE.Line( geometry, lineMaterial );
      line.type = THREE.LinePieces;
      toplevel.add( line );

      var optFlag = ( opt.indexOf('colz') != -1 || opt.indexOf('col') != -1 );
      var fcolor = d3.rgb(root_colors[histo['fFillColor']]);
      var fillcolor = new THREE.Color( 0xDDDDDD );
      fillcolor.setRGB(fcolor.r/255, fcolor.g/255, fcolor.b/255);
      var bin, wei;
      for ( i = 0; i < bins.length; ++i ) {
         wei = tz( optFlag ? maxbin : bins[i].z );
         bin = THREE.SceneUtils.createMultiMaterialObject(
            new THREE.CubeGeometry( 2*size/nbinsx, wei, 2*size/nbinsy ),
            [ new THREE.MeshLambertMaterial( { color: fillcolor.getHex(), shading: THREE.NoShading } ),
              wireMaterial ] );
         bin.position.x = tx( bins[i].x + (scalex/2));
         bin.position.y = wei/2;
         bin.position.z = -(ty( bins[i].y + (scaley/2)));
         bin.name = "x: [" + bins[i].x.toPrecision(4) + ", " + (bins[i].x + scalex).toPrecision(4) + "]<br>" +
                    "y: [" + bins[i].y.toPrecision(4) + ", " + (bins[i].y + scaley).toPrecision(4) + "]<br>" +
                    "entries: " + bins[i].z.toFixed();
         toplevel.add( bin );
      }
      // create a point light
      var pointLight = new THREE.PointLight( 0xcfcfcf );
      pointLight.position.set(0, 50, 250);
      scene.add(pointLight);

      //var directionalLight = new THREE.DirectionalLight( 0x7f7f7f );
      //directionalLight.position.set( 0, -70, 100 ).normalize();
      //scene.add( directionalLight );

      var camera = new THREE.PerspectiveCamera( 45, w / h, 1, 1000 );
      camera.position.set( 0, size/2, 500 );
      camera.lookat = cube;

      var renderer = Detector.webgl ? new THREE.WebGLRenderer( { antialias: true } ) :
                     new THREE.CanvasRenderer( { antialias: true } );
      renderer.setSize( w, h );
      $( vis[0][0] ).hide().parent().append( renderer.domElement );
      renderer.render( scene, camera );

      this.add3DInteraction(renderer, scene, camera, toplevel);
      return renderer;
   }

   JSROOTPainter.drawHistogram3D = function(vis, pad, histo, hframe) {
      var i, j, k, logx = false, logy = false, logz = false,
          gridx = false, gridy = false, gridz = false;
      var opt = histo['fOption'].toLowerCase();
      if (pad && typeof(pad) != 'undefined') {
         logx = pad['fLogx'];
         logy = pad['fLogy'];
         logz = pad['fLogz'];
         gridx = pad['fGridx'];
         gridy = pad['fGridy'];
         gridz = pad['fGridz'];
      }
      var fillcolor = root_colors[histo['fFillColor']];
      var linecolor = root_colors[histo['fLineColor']];
      if (histo['fFillColor'] == 0) {
         fillcolor = '#4572A7';
      }
      if (histo['fLineColor'] == 0) {
         linecolor = '#4572A7';
      }
      var nbinsx = histo['fXaxis']['fNbins'];
      var nbinsy = histo['fYaxis']['fNbins'];
      var nbinsz = histo['fZaxis']['fNbins'];
      var scalex = (histo['fXaxis']['fXmax'] - histo['fXaxis']['fXmin']) /
                    histo['fXaxis']['fNbins'];
      var scaley = (histo['fYaxis']['fXmax'] - histo['fYaxis']['fXmin']) /
                    histo['fYaxis']['fNbins'];
      var scalez = (histo['fZaxis']['fXmax'] - histo['fZaxis']['fXmin']) /
                    histo['fZaxis']['fNbins'];
      var maxbin = -1e32, minbin = 1e32;
      maxbin = d3.max(histo['fArray']);
      minbin = d3.min(histo['fArray']);
      var bins = new Array();
      for (i=0; i<=nbinsx+2; ++i) {
         for (var j=0; j<nbinsy+2; ++j) {
            for (var k=0; k<nbinsz+2; ++k) {
               var bin_content = histo.getBinContent(i, j, k);
               if (bin_content > minbin) {
                  var point = {
                     x:histo['fXaxis']['fXmin'] + (i*scalex),
                     y:histo['fYaxis']['fXmin'] + (j*scaley),
                     z:histo['fZaxis']['fXmin'] + (k*scalez),
                     n:bin_content
                  };
                  bins.push(point);
               }
            }
         }
      }
      var w = vis.attr("width"), h = vis.attr("height"), size = 100;
      if (logx) {
         var tx = d3.scale.log().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([-size, size]);
         var utx = d3.scale.log().domain([-size, size]).range([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]);
      } else {
         var tx = d3.scale.linear().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([-size, size]);
         var utx = d3.scale.linear().domain([-size, size]).range([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]);
      }
      if (logy) {
         var ty = d3.scale.log().domain([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]).range([-size, size]);
         var uty = d3.scale.log().domain([size, -size]).range([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]);
      } else {
         var ty = d3.scale.linear().domain([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]).range([-size, size]);
         var uty = d3.scale.linear().domain([size, -size]).range([histo['fYaxis']['fXmin'], histo['fYaxis']['fXmax']]);
      }
      if (logz) {
         var tz = d3.scale.log().domain([histo['fZaxis']['fXmin'], histo['fZaxis']['fXmax']]).range([-size, size]);
         var utz = d3.scale.log().domain([-size, size]).range([histo['fZaxis']['fXmin'], histo['fZaxis']['fXmax']]);
      } else {
         var tz = d3.scale.linear().domain([histo['fZaxis']['fXmin'], histo['fZaxis']['fXmax']]).range([-size, size]);
         var utz = d3.scale.linear().domain([-size, size]).range([histo['fZaxis']['fXmin'], histo['fZaxis']['fXmax']]);
      }

      // three.js 3D drawing
      var scene = new THREE.Scene();

      var toplevel = new THREE.Object3D();
      toplevel.rotation.x = 30 * Math.PI / 180;
      toplevel.rotation.y = 30 * Math.PI / 180;
      scene.add( toplevel );

      var wireMaterial = new THREE.MeshBasicMaterial( {
         color: 0x000000,
         wireframe: true,
         wireframeLinewidth: 0.5,
         side: THREE.DoubleSide } );

      // create a new mesh with cube geometry
      var cube = new THREE.Mesh( new THREE.CubeGeometry( size*2, size*2, size*2 ), wireMaterial);

      // add the cube to the scene
      toplevel.add( cube );

      var textMaterial = new THREE.MeshBasicMaterial( { color: 0x000000 } );

      // add the calibration vectors and texts
      var geometry = new THREE.Geometry();
      var imax, istep, len = 3, plen, sin45 = Math.sin(45);
      var text3d, text;
      var xmajors = tx.ticks(5);
      var xminors = tx.ticks(25);
      for ( i=-size, j=0, k=0; i<=size; ++i ) {
         var is_major = ( utx( i ) <= xmajors[j] && utx( i+1 ) > xmajors[j] ) ? true : false;
         var is_minor = ( utx( i ) <= xminors[k] && utx( i+1 ) > xminors[k] ) ? true : false;
         plen = ( is_major ? len + 2 : len) * sin45;
         if ( is_major ) {
            text3d = new THREE.TextGeometry( xmajors[j], {
               size: 7,
               height: 0,
               curveSegments: 10
            });
            ++j;

            text3d.computeBoundingBox();
            var centerOffset = 0.5 * ( text3d.boundingBox.max.x - text3d.boundingBox.min.x );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( i-centerOffset, -size-13, size+plen );
            toplevel.add( text );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( i+centerOffset, -size-13, -size-plen );
            text.rotation.y = Math.PI;
            toplevel.add( text );
         }
         if ( is_major || is_minor ) {
            ++k;
            geometry.vertices.push( new THREE.Vector3( i, -size, size ) );
            geometry.vertices.push( new THREE.Vector3( i, -size-plen, size+plen ) );
            geometry.vertices.push( new THREE.Vector3( i, -size, -size ) );
            geometry.vertices.push( new THREE.Vector3( i, -size-plen, -size-plen ) );
         }
      }
      var ymajors = ty.ticks(5);
      var yminors = ty.ticks(25);
      for ( i=size, j=0, k=0; i>-size; --i ) {
         var is_major = ( uty( i ) <= ymajors[j] && uty( i-1 ) > ymajors[j] ) ? true : false;
         var is_minor = ( uty( i ) <= yminors[k] && uty( i-1 ) > yminors[k] ) ? true : false;
         plen = ( is_major ? len + 2 : len) * sin45;
         if ( is_major ) {
            text3d = new THREE.TextGeometry( ymajors[j], {
               size: 7,
               height: 0,
               curveSegments: 10
            });
            ++j;

            text3d.computeBoundingBox();
            var centerOffset = 0.5 * ( text3d.boundingBox.max.x - text3d.boundingBox.min.x );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( size+plen, -size-13, i+centerOffset );
            text.rotation.y = Math.PI/2;
            toplevel.add( text );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( -size-plen, -size-13, i-centerOffset );
            text.rotation.y = -Math.PI/2;
            toplevel.add( text );
         }
         if ( is_major || is_minor ) {
            ++k;
            geometry.vertices.push( new THREE.Vector3( size, -size, i ) );
            geometry.vertices.push( new THREE.Vector3( size+plen, -size-plen, i ) );
            geometry.vertices.push( new THREE.Vector3( -size, -size, i ) );
            geometry.vertices.push( new THREE.Vector3( -size-plen, -size-plen, i ) );
         }
      }
      var zmajors = tz.ticks(5);
      var zminors = tz.ticks(25);
      for ( i=-size, j=0, k=0; i<=size; ++i ) {
         var is_major = ( utz( i ) <= zmajors[j] && utz( i+1 ) > zmajors[j] ) ? true : false;
         var is_minor = ( utz( i ) <= zminors[k] && utz( i+1 ) > zminors[k] ) ? true : false;
         plen = ( is_major ? len + 2 : len) * sin45;
         if ( is_major ) {
            text3d = new THREE.TextGeometry( zmajors[j], {
               size: 7,
               height: 0,
               curveSegments: 10
            });
            ++j;

            text3d.computeBoundingBox();
            var offset = 0.6 * ( text3d.boundingBox.max.x - text3d.boundingBox.min.x );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( size+offset+7, i-2.5, size+offset+7 );
            text.rotation.y = Math.PI*3/4;
            toplevel.add( text );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( size+offset+7, i-2.5, -size-offset-7 );
            text.rotation.y = -Math.PI*3/4;
            toplevel.add( text );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( -size-offset-7, i-2.5, size+offset+7 );
            text.rotation.y = Math.PI/4;
            toplevel.add( text );

            text = new THREE.Mesh( text3d, textMaterial );
            text.position.set( -size-offset-7, i-2.5, -size-offset-7 );
            text.rotation.y = -Math.PI/4;
            toplevel.add( text );
         }
         if ( is_major || is_minor ) {
            ++k;
            geometry.vertices.push( new THREE.Vector3( size, i, size ) );
            geometry.vertices.push( new THREE.Vector3( size+plen, i, size+plen ) );
            geometry.vertices.push( new THREE.Vector3( size, i, -size ) );
            geometry.vertices.push( new THREE.Vector3( size+plen, i, -size-plen ) );
            geometry.vertices.push( new THREE.Vector3( -size, i, size ) );
            geometry.vertices.push( new THREE.Vector3( -size-plen, i, size+plen ) );
            geometry.vertices.push( new THREE.Vector3( -size, i, -size ) );
            geometry.vertices.push( new THREE.Vector3( -size-plen, i, -size-plen ) );
         }
      }

      // add the calibration lines
      var lineMaterial = new THREE.LineBasicMaterial( { color: 0x000000 } );
      var line = new THREE.Line( geometry, lineMaterial );
      line.type = THREE.LinePieces;
      toplevel.add( line );

      // create the bin cubes
      var constx = (size*2 / histo['fXaxis']['fNbins']) / maxbin;
      var consty = (size*2 / histo['fYaxis']['fNbins']) / maxbin;
      var constz = (size*2 / histo['fZaxis']['fNbins']) / maxbin;

      var optFlag = ( opt.indexOf('colz') != -1 || opt.indexOf('col') != -1 );
      var fcolor = d3.rgb(root_colors[histo['fFillColor']]);
      var fillcolor = new THREE.Color( 0xDDDDDD );
      fillcolor.setRGB(fcolor.r/255, fcolor.g/255, fcolor.b/255);
      var bin, wei;
      for ( i = 0; i < bins.length; ++i ) {
         wei = ( optFlag ? maxbin : bins[i].n );
         if (opt.indexOf('box1') != -1) {
            bin = new THREE.Mesh( new THREE.SphereGeometry( 0.5 * wei * constx /*, 16, 16 */ ),
                      new THREE.MeshPhongMaterial( { color: fillcolor.getHex(),
                              specular: 0xbfbfbf /*, shading: THREE.NoShading */ } ) );
         }
         else {
            bin = THREE.SceneUtils.createMultiMaterialObject(
               new THREE.CubeGeometry( wei * constx, wei * constz, wei * consty ),
               [ new THREE.MeshLambertMaterial( { color: fillcolor.getHex(), shading: THREE.NoShading } ),
                 wireMaterial ] );
         }
         bin.position.x = tx( bins[i].x - (scalex/2));
         bin.position.y = tz( bins[i].z - (scalez/2));
         bin.position.z = -(ty( bins[i].y - (scaley/2)));
         bin.name = "x: [" + bins[i].x.toPrecision(4) + ", " + (bins[i].x + scalex).toPrecision(4) + "]<br>" +
                    "y: [" + bins[i].y.toPrecision(4) + ", " + (bins[i].y + scaley).toPrecision(4) + "]<br>" +
                    "z: [" + bins[i].z.toPrecision(4) + ", " + (bins[i].z + scalez).toPrecision(4) + "]<br>" +
                    "entries: " + bins[i].n.toFixed();
         toplevel.add( bin );
      }
      // create a point light
      var pointLight = new THREE.PointLight( 0xcfcfcf );
      pointLight.position.set(0, 50, 250);
      scene.add(pointLight);

      //var directionalLight = new THREE.DirectionalLight( 0x7f7f7f );
      //directionalLight.position.set( 0, -70, 100 ).normalize();
      //scene.add( directionalLight );

      var camera = new THREE.PerspectiveCamera( 45, w / h, 1, 1000 );
      camera.position.set( 0, 0, 500 );
      camera.lookat = cube;

      var renderer = Detector.webgl ? new THREE.WebGLRenderer( { antialias: true } ) :
                     new THREE.CanvasRenderer( { antialias: true } );
      renderer.setSize( w, h );
      $( vis[0][0] ).hide().parent().append( renderer.domElement );
      renderer.render( scene, camera );

      this.add3DInteraction(renderer, scene, camera, toplevel);
   }

   JSROOTPainter.drawHStack = function(vis, pad, stack, hframe) {
      // paint the list of histograms
      // By default, histograms are shown stacked.
      //    -the first histogram is paint
      //    -then the sum of the first and second, etc
      if (!'fHists' in stack) return;
      if (stack['fHists'].length == 0) return;
      var histos = stack['fHists'];
      var nhists = stack['fHists'].length;
      var opt = "";
      if ('fOption' in stack) opt = stack['fOption'].toLowerCase();
      var lsame = false;
      if (opt.indexOf("same") != -1) {
         lsame = true;
         opt.replace("same", "");
      }
      // compute the min/max of each axis
      var i, h;
      var xmin = 1e100;
      var xmax = -xmin;
      var ymin = 1e100;
      var ymax = -xmin;
      for (var i=0; i<nhists; ++i) {
         h = histos[i];
         if (h['fXaxis']['fXmin'] < xmin) xmin = h['fXaxis']['fXmin'];
         if (h['fXaxis']['fXmax'] > xmax) xmax = h['fXaxis']['fXmax'];
         if (h['fYaxis']['fXmin'] < ymin) ymin = h['fYaxis']['fXmin'];
         if (h['fYaxis']['fXmax'] > ymax) ymax = h['fYaxis']['fXmax'];
      }
      var nostack = opt.indexOf("nostack") == -1 ? false : true;
      if (opt.indexOf("nostack") == -1) stack.buildStack();
      var themin, themax;
      if (stack['fMaximum'] == -1111) themax = stack.getMaximum(opt);
      else themax = stack['fMaximum'];
      if (stack['fMinimum'] == -1111) {
         themin = stack.getMinimum(opt);
         if (pad && pad['fLogy']) {
            if (themin > 0) themin *= .9;
            else themin = themax * 1.e-3;
         }
         else if (themin > 0)
            themin = 0;
      }
      else themin = stack['fMinimum'];
      if (!('fHistogram' in stack)) {
         h = stack['fHists'][0];
         stack['fHistogram'] = new Object();
         stack['fHistogram']['_typename'] = stack['fHists'][0]['_typename'];
         JSROOTCore.addMethods(stack['fHistogram']);
         stack['fHistogram']['fName'] = "unnamed";
         stack['fHistogram']['fBits'] = 0;
         stack['fHistogram']['TestBit'] = function(bit) { return false };
         stack['fHistogram']['fOption'] = "";
         stack['fHistogram']['fXaxis'] = JSROOTCore.clone(h['fXaxis']);
         stack['fHistogram']['fYaxis'] = JSROOTCore.clone(h['fYaxis']);
         stack['fHistogram']['fXaxis']['fXmin'] = xmin;
         stack['fHistogram']['fXaxis']['fXmax'] = xmax;
         stack['fHistogram']['fYaxis']['fXmin'] = ymin;
         stack['fHistogram']['fYaxis']['fXmax'] = ymax;
         stack['fHistogram']['fXaxis']['fNbins'] = 0;
         stack['fHistogram']['fYaxis']['fNbins'] = 0;
         stack['fHistogram']['fArray'] = new Array();
      }
      stack['fHistogram']['fTitle'] = stack['fTitle'];
      //var histo = JSROOTCore.clone(stack['fHistogram']);
      var histo = stack['fHistogram'];
      if (!histo.TestBit(TH1StatusBits.kIsZoomed)) {
         if (nostack && stack['fMaximum'] != -1111) histo['fMaximum'] = stack['fMaximum'];
         else {
            if (pad && pad['fLogy']) histo['fMaximum'] = themax*(1+0.2*JSROOTMath.log10(themax/themin));
            else histo['fMaximum'] = 1.05 * themax;
         }
         if (nostack && stack['fMinimum'] != -1111) histo['fMinimum'] = stack['fMinimum'];
         else {
            if (pad && pad['fLogy']) histo['fMinimum'] = themin/(1+0.5*JSROOTMath.log10(themax/themin));
            else histo['fMinimum'] = themin;
         }
      }
      if (!lsame) {
         if (hframe) frame = hframe['frame'];
         else hframe = this.createFrame(vis, pad, histo, null);
         if ('fOption' in stack) {
            if (histo['fOption'].indexOf(stack['fOption'] == -1))
               histo['fOption'] += stack['fOption'];
         }
         if (histo['_typename'].match(/\bJSROOTIO.TH1/))
            this.drawHistogram1D(vis, pad, histo, hframe);
         else if (histo['_typename'].match(/\bJSROOTIO.TH2/))
            this.drawHistogram2D(vis, pad, histo, hframe);
         hframe['xmin'] = histo['fXaxis']['fXmin'];
         hframe['xmax'] = histo['fXaxis']['fXmax'];
         hframe['ymin'] = histo['fYaxis']['fXmin'];
         hframe['ymax'] = histo['fYaxis']['fXmax'];
      }
      if (nostack) {
         for (var i=0; i<nhists; ++i) {
            if ('redraw' in histo) // draw a clone if already drawn
               h = JSROOTCore.clone(histos[i]);
            else
               h = histos[i];
            if ('fOption' in stack) {
               if (h['fOption'].indexOf(stack['fOption'] == -1))
                  h['fOption'] += stack['fOption'];
            }
            h['fName'] += i;
            h['fOption'] += "same";
            // only draw TH1s (for the time being)
            if (h['_typename'].match(/\bJSROOTIO.TH1/))
               this.drawHistogram1D(vis, pad, h, hframe);
         }
      } else {
         var h1;
         for (var i=0; i<nhists; ++i) {
            if ('redraw' in histo) // draw a clone if already drawn
               h1 = JSROOTCore.clone(stack['fStack'][nhists-i-1]);
            else
               h1 = stack['fStack'][nhists-i-1];
            if ('fOption' in stack) {
               if (h1['fOption'].indexOf(stack['fOption'] == -1))
                  h1['fOption'] += stack['fOption'];
            }
            h1['fName'] += i;
            h1['fOption'] += "same";
            // only draw TH1s (for the time being)
            if (h1['_typename'].match(/\bJSROOTIO.TH1/))
               this.drawHistogram1D(vis, pad, h1, hframe);
         }
      }
   };

   JSROOTPainter.drawLatex = function(vis, string, x, y, attr) {
      var w = vis.attr("width"), h = vis.attr("height");
      while (string.indexOf('#') != -1)
         string = string.replace('#', '\\');
      string = string.replace(' ', '\\: ');

      var parse = new jsMath.Parser(string, null, null, null);
      parse.Parse();
      if (parse.error) return false;

      // method using jsMath do display formulae and LateX
      // unfortunately it works only on FireFox (Chrome displays it,
      // but at wrong coordinates, and IE doesn't support foreignObject
      // in SVG...)
      string = '\\displaystyle \\rm ' + string;
      var fo = vis.append("foreignObject")
         .attr("x", x)
         .attr("y", y)
         .attr("width", w - x)
         .attr("height", h - y);
      var math = fo.append("xhtml:div")
         .style("display", "inline")
         .style("color", attr['font-color'])
         .style('font-size', (attr['font-size']*0.98)+'px')
         .attr("class", "math")
         .html(string);
      jsMath.ProcessElement(math[0][0]);
      return true;
   };

   JSROOTPainter.drawLegend = function(vis, pad, pave) {
      var x=0, y=0, w=0, h=0;
      if (pave['fInit'] == 0) {
          x = pave['fX1'] * vis.attr("width")
          y = vis.attr("height") - pave['fY1'] * vis.attr("height");
          w = (pave['fX2'] - pave['fX1']) * vis.attr("width");
          h = (pave['fY2'] - pave['fY1']) * vis.attr("height");
      }
      else {
          x = pave['fX1NDC'] * vis.attr("width")
          y = vis.attr("height") - pave['fY1NDC'] * vis.attr("height");
          w = (pave['fX2NDC'] - pave['fX1NDC']) * vis.attr("width");
          h = (pave['fY2NDC'] - pave['fY1NDC']) * vis.attr("height");
      }
      y -= h;
      var fillcolor = root_colors[pave['fFillColor']];
      var lcolor = root_colors[pave['fLineColor']];
      var lwidth = pave['fBorderSize'] ? pave['fBorderSize'] : 0;
      if (pave['fFillStyle'] > 4000 && pave['fFillStyle'] < 4100)
         fillcolor = 'none';

      var p = vis.append("svg:g")
         .attr("width", w)
         .attr("height", h)
         .attr("transform", "translate(" + x + "," + y + ")");

      p.append("svg:rect")
         .attr("class", p)
         .attr("x", 0)
         .attr("y", 0)
         .attr("width", w)
         .attr("height", h)
         .attr("fill", fillcolor)
         .style("stroke-width", lwidth ? 1 : 0)
         .style("stroke", lcolor);

      var tcolor = root_colors[pave['fTextColor']];
      var tpos_x = pave['fMargin'] * w;
      var nlines = pave['fPrimitives'].length;
      var font_size = Math.round(h / (nlines * 1.5));
      //var font_size = Math.round(pave['fTextSize'] * vis.height());
      var fontDetails = getFontDetails(root_fonts[Math.floor(pave['fTextFont']/10)]);

      var max_len = 0, mul = 1.4;
      for (var j=0; j<nlines; ++j) {
         line = this.translateLaTeX(pave['fPrimitives'][j]['fLabel']);
         lw = tpos_x + stringWidth(vis, line, fontDetails['name'], fontDetails['weight'],
                                   font_size, fontDetails['style']);
         if (lw > max_len) max_len = lw;
      }
      if (max_len > w) {
         font_size *= 0.85 * (w / max_len);
         mul *= 0.95 * (max_len / w);
      }
      var x1 = pave['fX1NDC'];
      var x2 = pave['fX2NDC'];
      var y1 = pave['fY1NDC'];
      var y2 = pave['fY2NDC'];
      var margin = pave['fMargin']*( x2-x1 )/pave['fNColumns'];
      var yspace = (y2-y1)/nlines;
      var ytext = y2 + 0.5*yspace;  // y-location of 0th entry
      var boxw = margin*0.35;

      for (var i=0; i<nlines; ++i) {
         var leg = pave['fPrimitives'][i];
         var string = leg['fLabel'];
         var pos_y = ((i+1) * (font_size * mul)) - (font_size/3);
         var tpos_y = (i+1) * (font_size * mul);
         if (nlines == 1) {
            var pos_y = (h * 0.75) - (font_size/3);
            var tpos_y = h * 0.75;
         }

         var mo = gFile.GetMappedObject(leg['fObject']);
         if (mo) {
            leg['fFillColor']   = mo['fFillColor'];
            leg['fFillStyle']   = mo['fFillStyle'];
            leg['fLineColor']   = mo['fLineColor'];
            leg['fLineStyle']   = mo['fLineStyle'];
            leg['fLineWidth']   = mo['fLineWidth'];
            leg['fMarkerColor'] = mo['fMarkerColor'];
            leg['fMarkerSize']  = mo['fMarkerSize'];
            leg['fMarkerStyle'] = mo['fMarkerStyle'];
         }
         var line_color = root_colors[leg['fLineColor']];
         var line_width = leg['fLineWidth'];
         var line_style = root_line_styles[leg['fLineStyle']];
         var fill_color = root_colors[leg['fFillColor']];
         var fill_style = leg['fFillStyle'];
         var opt = leg['fOption'].toLowerCase();

         p.append("text")
            .attr("class", "text")
            .attr("text-anchor", "start")
            .attr("x", tpos_x)
            .attr("y", tpos_y)
            .attr("font-weight", fontDetails['weight'])
            .attr("font-style", fontDetails['style'])
            .attr("font-family", fontDetails['name'])
            .attr("font-size", font_size)
            .attr("fill", tcolor)
            .text(string);

         // Draw fill pattern (in a box)
         if (opt.indexOf('f') != -1) {
            // box total height is yspace*0.7
            // define x,y as the center of the symbol for this entry
            var xsym = margin/2;
            var ysym = ytext;
            var xf = new Array(4), yf = new Array(4);
            xf[0] = xsym - boxw;
            yf[0] = ysym - yspace*0.35;
            xf[1] = xsym + boxw;
            yf[1] = yf[0];
            xf[2] = xf[1];
            yf[2] = ysym + yspace*0.35;
            xf[3] = xf[0];
            yf[3] = yf[2];
            for (var j=0;j<4;j++) {
               xf[j] = xf[j] * vis.attr("width");
               yf[j] = yf[j] * vis.attr("height");
            }
            var ww = xf[1] - xf[0];
            var hh = yf[2] - yf[0];
            pos_y = pos_y - (hh/2);
            var pos_x = (tpos_x/2) - (ww/2);

            if (fill_style > 3000) {
               createFillPatterns(vis, fill_style, leg['fFillColor']);
               p.append("svg:rect")
                  .attr("x", pos_x)
                  .attr("y", pos_y)
                  .attr("width", ww)
                  .attr("height", hh)
                  .style("fill", "url(#pat" + fill_style + "_" + leg['fFillColor'] + ")")
                  .style("stroke-width", line_width)
                  .style("stroke", line_color);
            }
            else {
               p.append("svg:rect")
                  .attr("x", pos_x)
                  .attr("y", pos_y)
                  .attr("width", ww)
                  .attr("height", hh)
                  .attr("fill", fill_color)
                  .style("stroke-width", line_width)
                  .style("stroke", line_color);
            }
         }
         // Draw line
         if (opt.indexOf('l') != -1) {
            // line total length (in x) is margin*0.8
            var line_length = (0.7 * pave['fMargin']) * w;
            var pos_x = (tpos_x - line_length) / 2;
            p.append("svg:line")
               .attr("x1", pos_x)
               .attr("y1", pos_y)
               .attr("x2", pos_x + line_length)
               .attr("y2", pos_y)
               .style("stroke", line_color)
               .style("stroke-width", line_width)
               .style("stroke-dasharray", line_style);

         }
         // Draw error only
         if (opt.indexOf('e') != -1 && (opt.indexOf('l') == -1 || opt.indexOf('f') != -1)) {
         }
         // Draw Polymarker
         if (opt.indexOf('p') != -1) {

            var line_length = (0.7 * pave['fMargin']) * w;
            var pos_x = tpos_x / 2;

            var filled = false;
            if ((leg['fMarkerStyle'] == 8) ||
                (leg['fMarkerStyle'] > 19 && leg['fMarkerStyle'] < 24) ||
                (leg['fMarkerStyle'] == 29))
               filled = true;

            var info_marker = getRootMarker(root_markers, leg['fMarkerStyle']);

            var shape = info_marker['shape'];
            var filled = info_marker['toFill'];
            var toRotate = info_marker['toRotate'];
            var markerSize = leg['fMarkerSize'];
            var markerScale = 65;
            if (leg['fMarkerStyle'] == 1) markerScale = 1;

            switch (shape) {
               case 6:
                  var marker = "M " + (-4 * markerSize) + " " + (-1 * markerSize)
                              + " L " + 4 * markerSize + " " + (-1 * markerSize)
                              + " L " + (-2.4 * markerSize) + " " + 4 * markerSize
                              + " L 0 " + (-4 * markerSize) + " L " + 2.8 * markerSize
                              + " " + 4 * markerSize + " z";
                  break;
               case 7:
                  var marker = "M " + (- 4 * markerSize) + " " + (-4 * markerSize)
                              + " L " + 4 * markerSize + " " + 4 * markerSize + " M 0 "
                              + (-4 * markerSize) + " 0 " + 4 * markerSize + " M "
                              + 4 * markerSize + " " + (-4 * markerSize) + " L "
                              + (-4 * markerSize) + " " + 4 * markerSize + " M "
                              + (-4 * markerSize) + " 0 L " + 4 * markerSize + " 0";
                  break;
               default:
                  var marker = d3.svg.symbol()
                              .type(d3.svg.symbolTypes[shape])
                              .size(markerSize * markerScale);
                  break;
            }
            p.append("svg:path")
               .attr("transform", function(d) { return "translate(" + pos_x + "," + pos_y + ")"; })
               .style("fill", filled ? root_colors[leg['fMarkerColor']] : "none")
               .style("stroke", root_colors[leg['fMarkerColor']])
               .attr("d", marker);
         }
      }
      if (lwidth && lwidth > 1) {
         p.append("svg:line")
            .attr("x1", w+(lwidth/2))
            .attr("y1", lwidth+1)
            .attr("x2", w+(lwidth/2))
            .attr("y2", h+lwidth-1)
            .style("stroke", lcolor)
            .style("stroke-width", lwidth);
         p.append("svg:line")
            .attr("x1", lwidth+1)
            .attr("y1", h+(lwidth/2))
            .attr("x2", w+lwidth-1)
            .attr("y2", h+(lwidth/2))
            .style("stroke", lcolor)
            .style("stroke-width", lwidth);
      }
      return p;
   };

   JSROOTPainter.drawMultiGraph = function(vis, pad, mgraph, hframe) {
      var i, maximum, minimum, rwxmin=0, rwxmax=0, rwymin=0, rwymax=0, uxmin=0, uxmax=0, dx, dy;
      var npt = 100;
      var histo = mgraph['fHistogram'];
      var graphs = mgraph['fGraphs'];
      var scalex = 1, scaley = 1;
      var logx = false, logy = false, logz = false, gridx = false, gridy = false;
      var draw_all = true;
      if (pad && typeof(pad) != 'undefined') {
         rwxmin = pad.fUxmin;
         rwxmax = pad.fUxmax;
         rwymin = pad.fUymin;
         rwymax = pad.fUymax;
         logx = pad['fLogx']; logy = pad['fLogy']; logz = pad['fLogz'];
         gridx = pad['fGridx']; gridy = pad['fGridy'];
      }
      if ('fHistogram' in mgraph && mgraph['fHistogram']) {
         minimum = mgraph['fHistogram']['fYaxis']['fXmin'];
         maximum = mgraph['fHistogram']['fYaxis']['fXmax'];
         if (pad && typeof(pad) != 'undefined') {
            uxmin   = this.padtoX(pad, rwxmin);
            uxmax   = this.padtoX(pad, rwxmax);
         }
      } else {
         var g = graphs[0];
         if (g) {
            var r = g.computeRange();
            rwxmin = r['xmin']; rwymin = r['ymin'];
            rwxmax = r['xmax']; rwymax = r['ymax'];
         }
         for (i=1; i<graphs.length; ++i) {
            var rx1,ry1,rx2,ry2;
            g = graphs[i];
            var r = g.computeRange();
            rx1 = r['xmin']; ry1 = r['ymin'];
            rx2 = r['xmax']; ry2 = r['ymax'];
            if (rx1 < rwxmin) rwxmin = rx1;
            if (ry1 < rwymin) rwymin = ry1;
            if (rx2 > rwxmax) rwxmax = rx2;
            if (ry2 > rwymax) rwymax = ry2;
            if (g['fNpoints'] > npt) npt = g['fNpoints'];
         }
         if (rwxmin == rwxmax) rwxmax += 1.;
         if (rwymin == rwymax) rwymax += 1.;
         dx = 0.05*(rwxmax-rwxmin);
         dy = 0.05*(rwymax-rwymin);
         uxmin = rwxmin - dx;
         uxmax = rwxmax + dx;
         if (logy) {
            if (rwymin <= 0) rwymin = 0.001*rwymax;
            minimum = rwymin/(1+0.5*JSROOTMath.log10(rwymax/rwymin));
            maximum = rwymax*(1+0.2*JSROOTMath.log10(rwymax/rwymin));
         } else {
            minimum  = rwymin - dy;
            maximum  = rwymax + dy;
         }
         if (minimum < 0 && rwymin >= 0) minimum = 0;
         if (maximum > 0 && rwymax <= 0) maximum = 0;
      }
      if (mgraph['fMinimum'] != -1111) rwymin = minimum = mgraph['fMinimum'];
      if (mgraph['fMaximum'] != -1111) rwymax = maximum = mgraph['fMaximum'];
      if (uxmin < 0 && rwxmin >= 0) {
         if (logx) uxmin = 0.9*rwxmin;
         //else                 uxmin = 0;
      }
      if (uxmax > 0 && rwxmax <= 0) {
         if (logx) uxmax = 1.1*rwxmax;
         //else                 uxmax = 0;
      }
      if (minimum < 0 && rwymin >= 0) {
         if(logy) minimum = 0.9*rwymin;
         //else                minimum = 0;
      }
      if (maximum > 0 && rwymax <= 0) {
         if(logy) maximum = 1.1*rwymax;
         //else                maximum = 0;
      }
      if (minimum <= 0 && logy) minimum = 0.001*maximum;
      if (uxmin <= 0 && logx) {
         if (uxmax > 1000) uxmin = 1;
         else              uxmin = 0.001*uxmax;
      }
      rwymin = minimum;
      rwymax = maximum;
      if ('fHistogram' in mgraph && mgraph['fHistogram']) {
         mgraph['fHistogram']['fYaxis']['fXmin'] = rwymin;
         mgraph['fHistogram']['fYaxis']['fXmax'] = rwymax;
      }
      var frame;
      if (hframe) frame = hframe['frame'];
      else {
         hframe = this.createFrame(vis, pad, histo, null);
         frame = hframe['frame'];
      }
      // Create a temporary histogram to draw the axis
      if ('fHistogram' in mgraph && mgraph['fHistogram']) {
         this.drawHistogram1D(vis, pad, histo, hframe);
         hframe['xmin'] = histo['fXaxis']['fXmin'];
         hframe['xmax'] = histo['fXaxis']['fXmax'];
         hframe['ymin'] = histo['fYaxis']['fXmin'];
         hframe['ymax'] = histo['fYaxis']['fXmax'];
      }
      else {
         // the graph is created with at least as many channels as there are points
         // to permit zooming on the full range
         var w = vis.attr("width"), h = vis.attr("height");
         var label_font_size = Math.round(0.035 * h);
         w = frame.attr("width"); h = frame.attr("height");
         rwxmin = uxmin;
         rwxmax = uxmax;
         hframe['xmin'] = rwxmin;
         hframe['xmax'] = rwxmax;
         hframe['ymin'] = rwymin;
         hframe['ymax'] = rwymax;

         if (logx)
            var x = d3.scale.log().domain([rwxmin, rwxmax]).range([0, w]);
         else
            var x = d3.scale.linear().domain([rwxmin, rwxmax]).range([0, w]);
         if (logy)
            var y = d3.scale.log().domain([rwymin, rwymax]).range([h, 0]);
         else
            var y = d3.scale.linear().domain([rwymin, rwymax]).range([h, 0]);

         /* X-axis */
         var x_axis = d3.svg.axis()
            .scale(x)
            .orient("bottom")
            .tickPadding(5)
            .tickSubdivide(3)
            .tickSize(-0.03 * h, -0.03 * h / 2, null)
            .tickFormat(function(d,i) { return parseFloat(d.toPrecision(12)); })
            .ticks(10);

         /* Y-axis minor ticks */
         var y_axis = d3.svg.axis()
            .scale(y)
            .orient("left")
            .tickPadding(5)
            .tickSubdivide(3)
            .tickSize(-0.03 * w, -0.03 * w / 2, null)
            .tickFormat(function(d,i) { return parseFloat(d.toPrecision(12)); })
            .ticks(8);

         var xax = frame.append("svg:g")
            .attr("class", "xaxis")
            .attr("transform", "translate(0," + h + ")")
            .call(x_axis);

         var yax = frame.append("svg:g")
            .attr("class", "yaxis")
            .call(y_axis);

         var font_size = Math.round(0.050 * vis.attr("height"));
         vis.append("text")
            .attr("class", "title")
            .attr("text-anchor", "middle")
            .attr("x", vis.attr("width")/2)
            .attr("y", 1 + font_size)
            .attr("font-family", "Arial")
            .attr("font-size", font_size)
            .text(mgraph['fTitle']);

         xax.selectAll("text").attr("font-size", label_font_size);
         yax.selectAll("text").attr("font-size", label_font_size);

         frame['x_axis']  = x_axis;
         frame['y_axis']  = y_axis;
         frame['x_fsize'] = label_font_size;
         frame['y_fsize'] = label_font_size;
         frame['x_font']  = {'weight' : "",'style' : "", 'name' : "arial" };
         frame['y_font']  = {'weight' : "",'style' : "", 'name' : "arial" };

         var histo = new Object();
         histo['fXaxis'] = new Object();
         histo['fYaxis'] = new Object();
         histo['x'] = x;
         histo['y'] = y;
         histo['x_min'] = rwxmin;
         histo['x_max'] = rwxmax;
         histo['y_min'] = rwymin;
         histo['y_max'] = rwymax;
         histo['fLineWidth'] = 1;
         histo['fXaxis']['fNdivisions'] = 510;
         histo['fYaxis']['fNdivisions'] = 510;
         histo['redraw'] = function() {
            JSROOTPainter.drawGrid(frame, histo, pad, x, y);
         };
         histo.redraw();
         this.addInteraction(frame, histo);
      }
      for (var i=0; i<graphs.length; ++i) {
         graphs[i]['fName'] += i;
         this.drawGraph(vis, pad, graphs[i], hframe);
      }
   };

   JSROOTPainter.drawObject = function(obj, idx) {
      var i, svg = null;
      function draw(init) {

         var render_to = '#histogram' + idx;
         if (typeof($(render_to)[0]) == 'undefined') {
            obj = null;
            return;
         }
         $(render_to).empty();

         for (i=0; i<func_list.length; ++i) {
            func_list[i]['isDrawn'] = false;
         }
         if (obj['_typename'].match(/\bTCanvas/)) {
            svg = JSROOTPainter.drawCanvas(obj, idx);
            if (init == true)
               window.setTimeout(function() { $(render_to)[0].scrollIntoView(); }, 50);
            return;
         }
         svg = JSROOTPainter.createCanvas($(render_to), idx);
         if (svg == null) return false;
         if (obj['_typename'].match(/\bJSROOTIO.TH1/)) {
            JSROOTPainter.drawHistogram1D(svg, null, obj, null);
         }
         else if (obj['_typename'].match(/\bJSROOTIO.TH2/)) {
            var renderer = 0;
            var vid = 'view3d_' + obj['fName'];
            $('<div><input type="checkbox" id='+vid+' /><label for='+vid+'>View in 3D</label></div>')
               .css('padding', '10px').css('position', 'absolute').insertBefore( svg[0][0] );
            $('#'+vid).click(function(e) {
               if ( $(this).prop('checked') ) {
                  renderer = JSROOTPainter.drawHistogram2D3D(svg, null, obj, null);
               } else {
                  $( svg[0][0] ).show().parent().find( renderer.domElement ).remove();
               }
            });
            JSROOTPainter.drawHistogram2D(svg, null, obj, null);
         }
         else if (obj['_typename'].match(/\bJSROOTIO.TH3/)) {
            JSROOTPainter.drawHistogram3D(svg, null, obj, null);
         }
         else if (obj['_typename'].match(/\bJSROOTIO.THStack/)) {
            JSROOTPainter.drawHStack(svg, null, obj, null)
         }
         else if (obj['_typename'].match(/\bJSROOTIO.TProfile/)) {
            JSROOTPainter.drawProfile(svg, null, obj, null);
         }
         else if (obj['_typename'] == 'JSROOTIO.TF1') {
            JSROOTPainter.drawFunction(svg, null, obj, null);
         }
         else if (obj['_typename'].match(/\bTGraph/) ||
                  obj['_typename'].match(/\bRooHist/) ||
                  obj['_typename'].match(/\RooCurve/)) {
            JSROOTPainter.drawGraph(svg, null, obj, null);
         }
         else if (obj['_typename'] == 'JSROOTIO.TMultiGraph') {
            JSROOTPainter.drawMultiGraph(svg, null, obj, null);
         }
         else if (typeof(drawUserObject) == 'function') {
            drawUserObject(obj, svg);
         }
         if (init == true)
            window.setTimeout(function() { $(render_to)[0].scrollIntoView(); }, 50);
      };
      //A better idom for binding with resize is to debounce
      var debounce = function(fn, timeout) {
         var timeoutID = -1;
         return function() {
            if (timeoutID > -1) {
               window.clearTimeout(timeoutID);
            }
            timeoutID = window.setTimeout(fn, timeout);
         }
      };
      var debounced_draw = debounce(function() { draw(false); }, 100);
      $(window).resize(debounced_draw);
      draw(true);
   };

   JSROOTPainter.drawPad = function(vis, pad) {
      var width = vis.attr("width"), height = vis.attr("height");
      var x = pad['fAbsXlowNDC'] * width;
      var y = height - pad['fAbsYlowNDC'] * height;
      var w = pad['fAbsWNDC'] * width;
      var h = pad['fAbsHNDC'] * height;
      y -= h;

      var fillcolor = root_colors[pad['fFillColor']];
      if (pad['fFillStyle'] > 4000 && pad['fFillStyle'] < 4100)
         fillcolor = 'none';

      var border_width = pad['fLineWidth'];
      var border_color = root_colors[pad['fLineColor']];
      if (pad['fBorderMode'] == 0) {
         border_width = 0;
         border_color = 'none';
      }

      var new_pad = vis.append("svg:g")
         .attr("width", w)
         .attr("height", h)
         .attr("transform", "translate(" + x + "," + y + ")");

      new_pad.append("svg:rect")
         .attr("class", new_pad)
         .attr("x", 0)
         .attr("y", 0)
         .attr("width", w)
         .attr("height", h)
         .attr("fill", fillcolor)
         .style("stroke-width", border_width)
         .style("stroke", border_color);

      this.drawPrimitives(new_pad, pad);
      return new_pad;
   };

   JSROOTPainter.drawPaletteAxis = function(vis, palette, minbin, maxbin) {
      var width = vis.attr("width");
      var height = vis.attr("height");

      var pos_x = palette['fX1NDC'] * width;
      var pos_y = height - palette['fY1NDC'] * height;
      var s_width = Math.abs(palette['fX2NDC'] - palette['fX1NDC']) * width;
      var s_height = Math.abs(palette['fY2NDC'] - palette['fY1NDC']) * height;
      pos_y -= s_height;

      /*
       * Draw palette pad
       */
      var pal = vis.append("svg:g")
          .attr("height", s_height)
          .attr("width", s_width)
          .attr("transform", "translate(" + pos_x + ", " + pos_y + ")");

      var axis = palette['fAxis'];

      /*
       * Draw the default palette
       */
      var rectHeight = s_height / default_palette.length;
      pal.selectAll("colorRect")
         .data(default_palette)
         .enter()
         .append("svg:rect")
         .attr("class", "colorRect")
         .attr("x", 0)
         .attr("y", function(d, i) { return s_height - (i + 1)* rectHeight;})
         .attr("width", s_width)
         .attr("height", rectHeight)
         .attr("fill", function(d) { return d;})
         .attr("stroke", function(d) {return d;});
      /*
       * Build and draw axes
       */
      var nbr1 = 8;//Math.max((ndiv % 10000) % 100, 1);
      var nbr2 = 0;//Math.max(Math.round((ndiv % 10000) / 100), 1);
      var nbr3 = 0;//Math.max(Math.round(ndiv / 10000), 1);

      var z = d3.scale.linear().clamp(true)
               .domain([minbin, maxbin])
               .range([s_height, 0])
               .nice();

      var axisOffset = axis['fLabelOffset'] * width;
      var tickSize = axis['fTickSize'] * width;
      var z_axis = d3.svg.axis()
                  .scale(z)
                  .orient("right")
                  .tickSubdivide(nbr2)
                  .tickPadding(axisOffset)
                  .tickSize(-tickSize, -tickSize/2, 0)
                  .ticks(nbr1);

       var zax = pal.append("svg:g")
               .attr("class", "zaxis")
               .attr("transform", "translate(" + s_width + ", 0)")
               .call(z_axis);

      var axisFontDetails = getFontDetails(root_fonts[Math.floor(axis['fLabelFont'] /10)]);
      var axisLabelFontSize = axis['fLabelSize'] * height;
      zax.selectAll("text")
         .attr("font-size", axisLabelFontSize)
         .attr("font-weight", axisFontDetails['weight'])
         .attr("font-style", axisFontDetails['style'])
         .attr("font-family", axisFontDetails['name'])
         .attr("fill", root_colors[axis['fLabelColor']]);

      /*
       * Add palette axis title
       */
      var title = axis['fTitle'];
      if (title != "" && typeof(axis['fTitleFont']) != 'undefined') {
         axisFontDetails = getFontDetails(root_fonts[Math.floor(axis['fTitleFont'] /10)]);
         var axisTitleFontSize = axis['fTitleSize'] * height;
         pal.append("text")
               .attr("class", "Z axis label")
               .attr("x", s_width + axisLabelFontSize)
               .attr("y", s_height)
               .attr("text-anchor", "end")
               .attr("font-family", axisFontDetails['name'])
               .attr("font-weight", axisFontDetails['weight'])
               .attr("font-style", axisFontDetails['style'])
               .attr("font-size", axisTitleFontSize )
               .text(title);
      }
   };

   JSROOTPainter.drawPaveLabel = function(vis, pavelabel) {
      var w = vis.attr("width"), h = vis.attr("height");
      var pos_x = pavelabel['fX1NDC'] * w;
      var pos_y = (1.0 - pavelabel['fY1NDC']) * h;
      var width = Math.abs(pavelabel['fX2NDC'] - pavelabel['fX1NDC']) * w;
      var height = Math.abs(pavelabel['fY2NDC'] - pavelabel['fY1NDC']) * h;
      pos_y -= height;
      var font_size = Math.round(height / 1.9);
      var fcolor = root_colors[pavelabel['fFillColor']];
      var lcolor = root_colors[pavelabel['fLineColor']];
      var tcolor = root_colors[pavelabel['fTextColor']];
      var scolor = root_colors[pavelabel['fShadowColor']];
      if (pavelabel['fFillStyle'] == 0) fcolor = 'none';
      // align = 10*HorizontalAlign + VerticalAlign
      // 1=left adjusted, 2=centered, 3=right adjusted
      // 1=bottom adjusted, 2=centered, 3=top adjusted
      var align = 'start', halign = Math.round(pavelabel['fTextAlign']/10);
      var baseline = 'bottom', valign = pavelabel['fTextAlign']%10;
      if (halign == 1) align = 'start';
      else if (halign == 2) align = 'middle';
      else if (halign == 3) align = 'end';
      if (valign == 1) baseline = 'bottom';
      else if (valign == 2) baseline = 'middle';
      else if (valign == 3) baseline = 'top';
      var lmargin = 0;
      switch (halign) {
         case 1:
            lmargin = pavelabel['fMargin'] * width;
            break;
         case 2:
            lmargin = width/2;
            break;
         case 3:
            lmargin = width - (pavelabel['fMargin'] * width);
            break;
      }
      var lwidth = pavelabel['fBorderSize'] ? pavelabel['fBorderSize'] : 0;
      var fontDetails = getFontDetails(root_fonts[Math.floor(pavelabel['fTextFont']/10)]);

      var pave = vis.append("svg:g")
         .attr("width", width)
         .attr("height", height)
         .attr("transform", "translate(" + pos_x + "," + pos_y + ")");

      pave.append("svg:rect")
         .attr("class", pave)
         .attr("x", 0)
         .attr("y", 0)
         .attr("width", width)
         .attr("height", height)
         .attr("fill", fcolor)
         .style("stroke-width", lwidth ? 1 : 0)
         .style("stroke", lcolor);

      var line = this.translateLaTeX(pavelabel['fLabel']);

      var lw = stringWidth(vis, line, fontDetails['name'], fontDetails['weight'],
                           font_size, fontDetails['style']);
      if (lw > width)
         font_size *= 0.98 * (width / lw);

      pave.append("text")
         .attr("class", "text")
         .attr("text-anchor", align)
         .attr("x", lmargin)
         .attr("y", (height/2) + (font_size/3))
         .attr("font-weight", fontDetails['weight'])
         .attr("font-style", fontDetails['style'])
         .attr("font-family", fontDetails['name'])
         .attr("font-size", font_size)
         .attr("fill", tcolor)
         .text(line);

      if (lwidth && lwidth > 1) {
         pave.append("svg:line")
            .attr("x1", width+(lwidth/2))
            .attr("y1", lwidth+1)
            .attr("x2", width+(lwidth/2))
            .attr("y2", height+lwidth-1)
            .style("stroke", lcolor)
            .style("stroke-width", lwidth);
         pave.append("svg:line")
            .attr("x1", lwidth+1)
            .attr("y1", height+(lwidth/2))
            .attr("x2", width+lwidth-1)
            .attr("y2", height+(lwidth/2))
            .style("stroke", lcolor)
            .style("stroke-width", lwidth);
      }
   };

   JSROOTPainter.drawPaveText = function(vis, pavetext) {
      var i, j, lw, w = vis.attr("width"), h = vis.attr("height");

      var pos_x = pavetext['fX1NDC'] * w;
      var pos_y = (1.0 - pavetext['fY1NDC']) * h;
      var width = Math.abs(pavetext['fX2NDC'] - pavetext['fX1NDC']) * w;
      var height = Math.abs(pavetext['fY2NDC'] - pavetext['fY1NDC']) * h;
      pos_y -= height;
      var line, nlines = pavetext['fLines'].length;
      var font_size = Math.round(height / (nlines * 1.5));
      var fcolor = root_colors[pavetext['fFillColor']];
      var lcolor = root_colors[pavetext['fLineColor']];
      var tcolor = root_colors[pavetext['fTextColor']];
      var scolor = root_colors[pavetext['fShadowColor']];
      if (pavetext['fFillStyle'] == 0) fcolor = 'none';
      // align = 10*HorizontalAlign + VerticalAlign
      // 1=left adjusted, 2=centered, 3=right adjusted
      // 1=bottom adjusted, 2=centered, 3=top adjusted
      // "middle", "start", "end"

      var align = 'start', halign = Math.round(pavetext['fTextAlign']/10);
      var baseline = 'bottom', valign = pavetext['fTextAlign']%10;
      if (halign == 1) align = 'start';  else
      if (halign == 2) align = 'middle';  else
      if (halign == 3) align = 'end';
      if (valign == 1) baseline = 'bottom'; else
      if (valign == 2) baseline = 'middle'; else
      if (valign == 3) baseline = 'top';

      var lmargin = 0;
      switch (halign) {
         case 1:
            lmargin = pavetext['fMargin'] * width;
            break;
         case 2:
            lmargin = width/2;
            break;
         case 3:
            lmargin = width - (pavetext['fMargin'] * width);
            break;
      }


      // for now ignore all align parameters, draw as is
      if (nlines>1) lmargin = pavetext['fMargin'] * width / 2;

      var fontDetails = getFontDetails(root_fonts[Math.floor(pavetext['fTextFont']/10)]);
      var lwidth = pavetext['fBorderSize'] ? pavetext['fBorderSize'] : 0;

      // $("#report").append("<br> draw pave");

      this.draw_g = vis.append("svg:g")
         .attr("x", pos_x)
         .attr("y", pos_y)
         .attr("width", width)
         .attr("height", height)
         .attr("transform", "translate(" + pos_x + "," + pos_y + ")");

      this.draw_g.append("svg:rect")
         .attr("x", 0)
         .attr("y", 0)
         .attr("height", height)
         .attr("width", width)
         .attr("fill", fcolor)
         .style("stroke-width", lwidth ? 1 : 0)
         .style("stroke", lcolor);

      var first_stat = 0;
      var num_cols = 0;
      var max_len = 0;

      for (var j=0; j<nlines; ++j) {
         var line = pavetext['fLines'][j]['fTitle'];
         line = JSROOTPainter.translateLaTeX(line);
         var parts = line.split("|");
         if (parts && (parts.length>1)) {
            if (first_stat == 0) first_stat = j;
            if (parts.length > num_cols) num_cols = parts.length;
         }
         var lw = lmargin + stringWidth(vis, line, fontDetails['name'], fontDetails['weight'],
                   font_size, fontDetails['style']);
         if (lw > max_len) max_len = lw;
      }

      if (max_len > width)
         font_size *= 0.98 * (width / max_len);

      if (nlines == 1) {
         line = JSROOTPainter.translateLaTeX(pavetext['fLines'][0]['fTitle']);

         this.draw_g.append("text")
            .attr("text-anchor", align)
            .attr("x", lmargin)
            .attr("y", (height/2) + (font_size/3))
            .attr("font-family", fontDetails['name'])
            .attr("font-weight", fontDetails['weight'])
            .attr("font-style", fontDetails['style'])
            .attr("font-size", font_size)
            .attr("fill", tcolor)
            .text(line);
      }
      else {

         for (j=0; j<nlines; ++j) {
            var jcolor = root_colors[pavetext['fLines'][j]['fTextColor']];
            if (pavetext['fLines'][j]['fTextColor'] == 0)  jcolor = tcolor;
            line = JSROOTPainter.translateLaTeX(pavetext['fLines'][j]['fTitle']);
            if (pavetext['_typename'] == 'JSROOTIO.TPaveStats') {
               var posy = (j == 0) ? (font_size * 1.2) :
                                     font_size * (0.05 + (j+1)*1.4);
               if ((first_stat > 0) && (j >= first_stat)) {
                  posy -= font_size*0.3; // dut to middle allignment
                  var parts = line.split("|");
                  for (var n=0;n<parts.length;n++)
                     this.draw_g.append("text")
                     .attr("text-anchor", "middle")
                     .attr("x",  width*(n+0.5)/num_cols)
                     .attr("y", posy + font_size*0.6)
                     .attr("font-family", fontDetails['name'])
                     .attr("font-weight", fontDetails['weight'])
                     .attr("font-style", fontDetails['style'])
                     .attr("font-size", font_size)
                     .attr("fill", jcolor)
                     .text(parts[n]);
               } else
               if ((j==0) || (line.indexOf('=')<0)) {
                  this.draw_g.append("text")
                     .attr("text-anchor", (j == 0) ? "middle" : "start")
                     .attr("x", ((j==0) ? width/2 : pavetext['fMargin'] * width))
                     .attr("y", posy)
                     .attr("font-family", fontDetails['name'])
                     .attr("font-weight", fontDetails['weight'])
                     .attr("font-style", fontDetails['style'])
                     .attr("font-size", font_size)
                     .attr("fill", jcolor)
                     .text(line);
               } else {
                  var parts = line.split("=");
                  for (var n=0;n<2;n++)
                     this.draw_g.append("text")
                     .attr("text-anchor", (n==0) ? "start" : "end")
                     .attr("x", (n==0) ? pavetext['fMargin'] * width : (1-pavetext['fMargin']) * width)
                     .attr("y", posy)
                     .attr("font-family", fontDetails['name'])
                     .attr("font-weight", fontDetails['weight'])
                     .attr("font-style", fontDetails['style'])
                     .attr("font-size", font_size)
                     .attr("fill", jcolor)
                     .text(parts[n]);
               }
            } else {
               this.draw_g.append("text")
                  .attr("text-anchor", "start")
                  .attr("x", lmargin)
                  .attr("y", (j+1) * (font_size * 1.4))
                  .attr("font-family", fontDetails['name'])
                  .attr("font-weight", fontDetails['weight'])
                  .attr("font-style", fontDetails['style'])
                  .attr("font-size", font_size)
                  .attr("fill", jcolor)
                  .text(line);
            }
         }
      }

      if (pavetext['fBorderSize'] && pavetext['_typename'] == 'JSROOTIO.TPaveStats') {
         this.draw_g.append("svg:line")
            .attr("class", "pavedraw")
            .attr("x1", 0)
            .attr("y1", font_size * 1.6)
            .attr("x2", width)
            .attr("y2", font_size * 1.6)
            .style("stroke", lcolor)
            .style("stroke-width", lwidth ? 1 : 'none');
      }

      if ((first_stat > 0) && (num_cols > 1)) {
         var yy = (1.4 * first_stat + 0.6) * font_size;

         this.draw_g.append("svg:line")
         .attr("x1", 0)
         .attr("y1", yy)
         .attr("x2", width)
         .attr("y2", yy)
         .style("stroke", lcolor)
         .style("stroke-width", lwidth ? 1 : 'none');

         for (var ncol = 0; ncol<num_cols-1; ncol++)
            this.draw_g.append("svg:line")
            .attr("x1", width/num_cols*(ncol+1))
            .attr("y1", yy)
            .attr("x2", width/num_cols*(ncol+1) )
            .attr("y2", height)
            .style("stroke", lcolor)
            .style("stroke-width", lwidth ? 1 : 'none');

      }

      if (lwidth && lwidth > 1) {
         this.draw_g.append("svg:line")
            .attr("x1", width+(lwidth/2))
            .attr("y1", lwidth+1)
            .attr("x2", width+(lwidth/2))
            .attr("y2", height+lwidth-1)
            .style("stroke", lcolor)
            .style("stroke-width", lwidth);
         this.draw_g.append("svg:line")
            .attr("x1", lwidth+1)
            .attr("y1", height+(lwidth/2))
            .attr("x2", width+lwidth-1)
            .attr("y2", height+(lwidth/2))
            .style("stroke", lcolor)
            .style("stroke-width", lwidth);
      }
   };

   JSROOTPainter.drawPrimitives = function(vis, pad) {
      var i, j, fframe = null, frame = null;
      var primitives = pad['fPrimitives'];
      for (i=0; i<primitives.length; ++i) {
         var classname = primitives[i]['_typename'];
         if (classname == 'JSROOTIO.TFrame') {
            fframe = frame = this.createFrame(vis, pad, null, primitives[i]);
         }
         if (classname == 'JSROOTIO.TPad') {
            this.drawPad(vis, primitives[i])
         }
         if (classname == 'JSROOTIO.TPaveLabel') {
            this.drawPaveLabel(vis, primitives[i]);
         }
         if (classname == 'JSROOTIO.TLegend') {
            this.drawLegend(vis, pad, primitives[i]);
         }
         if (classname == 'JSROOTIO.TPaveText') {
            this.drawPaveText(vis, primitives[i]);
         }
         if (classname == 'JSROOTIO.TLatex') {
            this.drawText(vis, pad, primitives[i]);
         }
         if (classname == 'JSROOTIO.TText') {
            this.drawText(vis, pad, primitives[i]);
         }
         if (classname.match(/\bJSROOTIO.TH1/)) {
            this.drawHistogram1D(vis, pad, primitives[i], frame);
            if (fframe) {
               fframe['xmin'] = primitives[i]['fXaxis']['fXmin'];
               fframe['xmax'] = primitives[i]['fXaxis']['fXmax'];
               fframe['ymin'] = primitives[i]['fYaxis']['fXmin'];
               fframe['ymax'] = primitives[i]['fYaxis']['fXmax'];
               if (primitives[i]['fXaxis'].TestBit(EAxisBits.kAxisRange)) {
                  fframe['xmin'] = primitives[i].getBinLowEdge(primitives[i]['fXaxis']['fFirst']);
                  fframe['xmax'] = primitives[i].getBinUpEdge(primitives[i]['fXaxis']['fLast']);
               }
            }
         }
         if (classname.match(/\bJSROOTIO.TH2/)) {
            this.drawHistogram2D(vis, pad, primitives[i], frame);
         }
         if (classname.match(/\bJSROOTIO.THStack/)) {
            this.drawHStack(vis, pad, primitives[i], frame)
         }
         if (classname.match(/\bJSROOTIO.TProfile/)) {
            this.drawProfile(vis, pad, primitives[i], frame);
         }
         if (classname == 'JSROOTIO.TF1') {
            if (typeof(primitives[i]['isDrawn']) == 'undefined' || primitives[i]['isDrawn'] == false)
               this.drawFunction(vis, pad, primitives[i], fframe);
            primitives[i]['isDrawn'] = true;
         }
         if (classname.match(/\bTGraph/) ||
             classname.match(/\bRooHist/) ||
             classname.match(/\RooCurve/)) {
            primitives[i]['fName'] += i;
            this.drawGraph(vis, pad, primitives[i], frame);
         }
         if (classname == 'JSROOTIO.TMultiGraph') {
            this.drawMultiGraph(vis, pad, primitives[i], frame);
         }
      }
   };

   JSROOTPainter.drawProfile = function(vis, pad, histo, hframe) {
      var i, logx = false, logy = false, logz = false, gridx = false, gridy = false;
      if (pad && typeof(pad) != 'undefined') {
         logx = pad['fLogx'];
         logy = pad['fLogy'];
         logz = pad['fLogz'];
         gridx = pad['fGridx'];
         gridy = pad['fGridy'];
      }
      var fillcolor = root_colors[histo['fFillColor']];
      var linecolor = root_colors[histo['fLineColor']];
      if (histo['fFillColor'] == 0) {
         fillcolor = '#4572A7';
      }
      if (histo['fLineColor'] == 0) {
         linecolor = '#4572A7';
      }
      //histo['fgApproximate'] = true;
      var binwidth = ((histo['fXaxis']['fXmax'] - histo['fXaxis']['fXmin']) / histo['fXaxis']['fNbins']);
      var bins = d3.range(histo['fXaxis']['fNbins']).map(function(p) {
         return {
            x:  histo['fXaxis']['fXmin'] + (p * binwidth) + (binwidth / 2.0),
            y:  histo.getBinContent(p+1),
            xerr: binwidth / 2.0,
            yerr: histo.getBinError(p+1)
         };
      });
      var ret = hframe != null ? hframe : this.createFrame(vis, pad, histo, null);
      var frame = ret['frame'];
      var svg_frame = d3.select(ret['id']);
      var w = frame.attr("width"), h = frame.attr("height");
      if (logx)
         var x = d3.scale.log().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([0, w]);
      else
         var x = d3.scale.linear().domain([histo['fXaxis']['fXmin'], histo['fXaxis']['fXmax']]).range([0, w]);
      if (logy)
         var y = d3.scale.log().domain([histo['fYmin'], histo['fYmax']]).range([h, 0]);
      else
         var y = d3.scale.linear().domain([histo['fYmin'], histo['fYmax']]).range([h, 0]);

      histo['x_min'] = histo['fXaxis']['fXmin'];
      histo['x_max'] = histo['fXaxis']['fXmax'];
      histo['y_min'] = histo['fYmin'];
      histo['y_max'] = histo['fYmax'];

      histo['x'] = x;
      histo['y'] = y;
      histo['bins'] = bins;

      this.drawGrid(frame, histo, pad, x, y);
      this.drawErrors(svg_frame, bins, histo, pad, x, y);
      this.drawAxes(frame, histo, pad, x, y);
      this.drawTitle(vis, histo, pad);
      this.addInteraction(frame, histo);
      var draw_stats = this.drawFunctions(vis, histo, pad, ret);
      if (draw_stats && !pad || typeof(pad) == 'undefined')
         this.drawStat(vis, histo);
      return null;
   };

   JSROOTPainter.addLine = function(pavetext, txt) {
      pavetext['fLines'].push( {'fTitle': txt, "fTextColor": 1} );
   }

   JSROOTPainter.drawStat = function(vis, histo) {
      var stats = {};
      stats['_typename'] = 'JSROOTIO.TPaveStats';
      stats['fName'] = 'stats';
      
      stats['fX1NDC'] = 0.78;
      stats['fY1NDC'] = 0.75;
      stats['fX2NDC'] = 0.78 + 0.2;
      stats['fY2NDC'] = 0.75 + 0.16;
      if ((histo['_typename'] && histo['_typename'].match(/\bTProfile/)) ||
          (histo['_typename'] && histo['_typename'].match(/\bTH2/)))
         stats['fY1NDC'] = 0.67;

      stats['fOptStat'] = 1111;
      stats['fLongest'] = 17;
      stats['fMargin'] = 0.05;

      stats['fBorderSize'] = 1;
      stats['fInit'] = 1;
      stats['fShadowColor'] = 1;
      stats['fCornerRadius'] = 0;

      stats['fX1'] = 1;
      stats['fY1'] = 100;
      stats['fX2'] = 1;
      stats['fY2'] = 100;

      stats['fResizing'] = false;
      stats['fBits'] = 0x03000009;
      stats['fLineColor'] = 1;
      stats['fLineStyle'] = 1;
      stats['fLineWidth'] = 1;

      stats['fFillColor'] = 0;
      stats['fFillStyle'] = 1001;
      
      stats['fTextAngle'] = 0;
      stats['fTextSize'] = 0;
      stats['fTextAlign'] = 12;
      stats['fTextColor'] = 1;
      stats['fTextFont'] = 42;
      
      stats['fLines'] = new Array;
      stats['fLines'].push({'fTitle': "hname", "fTextColor": 1});
      stats['fLines'][0]['fTitle'] = histo['fName']; 

      this.addLine(stats, "Entries = " + histo['fEntries']);
      if (histo['_typename'] && histo['_typename'].match(/\bTH2/))
         this.addLine(stats, "Mean x = " + histo.getMean(1).toFixed(6.4));
      else
         this.addLine(stats, "Mean = " + histo.getMean(1).toFixed(6.4));
      if ((histo['_typename'] && histo['_typename'].match(/\bTProfile/)) ||
          (histo['_typename'] && histo['_typename'].match(/\bTH2/))) {
         this.addLine(stats, "Mean y = " + histo.getMean(2).toFixed(6.4));
      }
      if (histo['_typename'] && histo['_typename'].match(/\bTH2/))
         this.addLine(stats, "RMS x = " + histo.getRMS(1).toFixed(6.4));
      else
         this.addLine(stats, "RMS = " + histo.getRMS(1).toFixed(6.4));
      if ((histo['_typename'] && histo['_typename'].match(/\bTProfile/)) ||
          (histo['_typename'] && histo['_typename'].match(/\bTH2/))) {
         this.addLine(stats, "RMS y = " + histo.getRMS(2).toFixed(6.4));
      }
      this.drawPaveText(vis, stats);
   };

   JSROOTPainter.drawText = function(vis, pad, text) {
      // align = 10*HorizontalAlign + VerticalAlign
      // 1=left adjusted, 2=centered, 3=right adjusted
      // 1=bottom adjusted, 2=centered, 3=top adjusted
      var i, w = vis.attr("width"), h = vis.attr("height");
      var align = 'start', halign = Math.round(text['fTextAlign']/10);
      var baseline = 'bottom', valign = text['fTextAlign']%10;
      if (halign == 1) align = 'start';
      else if (halign == 2) align = 'middle';
      else if (halign == 3) align = 'end';
      if (valign == 1) baseline = 'bottom';
      else if (valign == 2) baseline = 'middle';
      else if (valign == 3) baseline = 'top';
      var lmargin = 0;
      switch (halign) {
         case 1:
            lmargin = text['fMargin'] * w;
            break;
         case 2:
            lmargin = w/2;
            break;
         case 3:
            lmargin = w - (text['fMargin'] * w);
            break;
      }
      var font_size = Math.round(text['fTextSize'] * 0.7 * h);
      if (text.TestBit(kTextNDC)) {
         var pos_x = pad['fX1'] + text['fX']*(pad['fX2'] - pad['fX1']);
         var pos_y = pad['fY1'] + text['fY']*(pad['fY2'] - pad['fY1']);
      }
      else {
         font_size = Math.round(text['fTextSize'] * h);
         var pos_x = this.xtoPad(text['fX'], pad);
         var pos_y = this.ytoPad(text['fY'], pad);
      }
      pos_x = ((Math.abs(pad['fX1'])+pos_x)/(pad['fX2'] - pad['fX1']))*w;
      pos_y = (1-((Math.abs(pad['fY1'])+pos_y)/(pad['fY2'] - pad['fY1'])))*h;
      var tcolor = root_colors[text['fTextColor']];
      var fontDetails = getFontDetails(root_fonts[Math.floor(text['fTextFont']/10)]);

      var string = text['fTitle'];
      // translate the LaTeX symbols
      if (text['_typename'] == 'JSROOTIO.TLatex')
         string = this.translateLaTeX(text['fTitle']);

      vis.append("text")
         .attr("class", "text")
         .attr("x", pos_x)
         .attr("y", pos_y)
         .attr("font-family", fontDetails['name'])
         .attr("font-weight", fontDetails['weight'])
         .attr("font-style", fontDetails['style'])
         .attr("font-size", font_size)
         .attr("text-anchor", align)
         .attr("fill", tcolor)
         .text(string);
   };

   JSROOTPainter.drawTitle = function(vis, histo, pad) {
      /* draw the title only if we don't draw from a pad (see Olivier for details) */
      var w = vis.attr("width"), h = vis.attr("height");
      var font_size = Math.round(0.050 * h);
      var l_title = this.translateLaTeX(histo['fTitle']);
      if (!pad || typeof(pad) == 'undefined') {
         vis.append("text")
            .attr("class", "title")
            .attr("text-anchor", "middle")
            .attr("x", w/2)
            .attr("y", 1 + font_size)
            .attr("font-family", "Arial")
            .attr("font-size", font_size)
            .text(l_title);
      }
   };

   /**
    * List tree (dtree) related functions
    */

   JSROOTPainter.displayBranches = function(branches, dir_id, k) {
      for (var i=0; i<branches.length; ++i) {
         var nb_leaves = branches[i]['fLeaves'].length;
         var disp_name = branches[i]['fName'];
         var node_img = source_dir+'img/branch.png';
         var node_title = disp_name;
         var tree_link = "";
         if (nb_leaves == 0) {
            node_img = source_dir+'img/leaf.png';
         }
         else if (nb_leaves == 1 && branches[i]['fLeaves'][0]['fName'] == disp_name) {
            node_img = source_dir+'img/leaf.png';
            nb_leaves--;
         }
         if (branches[i]['fBranches'].length > 0) {
            node_img = source_dir+'img/branch.png';
         }
         key_tree.add(k, dir_id, disp_name, tree_link, node_title, '', node_img, node_img);
         nid = k; k++;
         for (var j=0; j<nb_leaves; ++j) {
            var disp_name = branches[i]['fLeaves'][j]['fName'];
            var node_title = disp_name;
            var node_img = source_dir+'img/leaf.png';
            var tree_link = "";
            key_tree.add(k, nid, disp_name, tree_link, node_title, '', node_img, node_img);
            k++;
         }
         if (branches[i]['fBranches'].length > 0) {
            k = JSROOTPainter.displayBranches(branches[i]['fBranches'], nid, k);
         }
      }
      return k;
   };

   JSROOTPainter.displayTree = function(tree, container, dir_id) {
      var tree_link = '';
      var content = "<p><a href='javascript: key_tree.openAll();'>open all</a> | <a href='javascript: key_tree.closeAll();'>close all</a></p>";
      var k = key_tree.aNodes.length;
      JSROOTPainter.displayBranches(tree['fBranches'], dir_id, k);
      content += key_tree;
      $(container).append(content);
      key_tree.openTo(dir_id, true);
   };

   JSROOTPainter.displayListOfKeys = function(keys, container) {
      delete key_tree;
      var content = "<p><a href='javascript: key_tree.openAll();'>open all</a> | <a href='javascript: key_tree.closeAll();'>close all</a></p>";
      key_tree = new dTree('key_tree');
      key_tree.config.useCookies = false;
      key_tree.add(0, -1, 'File Content');
      var k = 1;
      var tree_link = '';
      for (var i=0; i<keys.length; ++i) {
         var message = keys[i]['className']+' is not yet implemented.';
         tree_link = "javascript:  alert('" + message + "')";
         var node_img = source_dir+'img/page.gif';
         if (keys[i]['className'].match(/\bTH1/) ||
             keys[i]['className'].match(/\bRooHist/)) {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/histo.png';
         }
         else if (keys[i]['className'].match(/\bTH2/)) {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/histo2d.png';
         }
         else if (keys[i]['className'].match(/\bTH3/)) {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/histo3d.png';
         }
         else if (keys[i]['className'].match(/\bTGraph/) ||
             keys[i]['className'].match(/\RooCurve/)) {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/graph.png';
         }
         else if (keys[i]['className'] ==  'TF1') {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/graph.png';
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'] ==  'TProfile') {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/profile.png';
            node_title = keys[i]['name'];
         }
         else if (keys[i]['name'] == 'StreamerInfo') {
            tree_link = "javascript: displayStreamerInfos(gFile.fStreamerInfo.fStreamerInfos);";
            node_img = source_dir+'img/question.gif';
         }
         else if (keys[i]['className'] == 'TDirectory' || keys[i]['className'] == 'TDirectoryFile') {
            tree_link = "javascript: showDirectory('"+keys[i]['name']+"',"+keys[i]['cycle']+","+(i+1)+");";
            node_img = source_dir+'img/folder.gif';
         }
         else if (keys[i]['className'] == 'TList' || keys[i]['className'] == 'TObjArray') {
            tree_link = "javascript: showCollection('"+keys[i]['name']+"',"+keys[i]['cycle']+","+(i+1)+");";
            node_img = source_dir+'img/folder.gif';
         }
         else if (keys[i]['className'] == 'TTree' || keys[i]['className'] == 'TNtuple') {
            tree_link = "javascript: readTree('"+keys[i]['name']+"',"+keys[i]['cycle']+","+(i+1)+");";
            node_img = source_dir+'img/tree.png';
         }
         else if (keys[i]['className'].match('TGeoManager') ||
                  keys[i]['className'].match('TGeometry')) {
            node_img = source_dir+'img/folder.gif';
         }
         else if (keys[i]['className'].match('TCanvas')) {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/canvas.png';
            node_title = keys[i]['name'];
         }
         if (keys[i]['name'] != '' && keys[i]['className'] != 'TFile')
            if (keys[i]['className'] == 'TDirectory' || keys[i]['className'] == 'TList' ||
                keys[i]['className'] == 'TObjArray' || keys[i]['className'] == 'TDirectoryFile')
               key_tree.add(k, 0, keys[i]['name']+';'+keys[i]['cycle'], tree_link, keys[i]['name'], '', node_img,
                            source_dir+'img/folderopen.gif');
            else if (keys[i]['className'] == 'TTree' || keys[i]['className'] == 'TNtuple')
               key_tree.add(k, 0, keys[i]['name']+';'+keys[i]['cycle'], tree_link, keys[i]['name'], '', node_img, node_img);
            else
               key_tree.add(k, 0, keys[i]['name']+';'+keys[i]['cycle'], tree_link, keys[i]['name'], '', node_img);
            k++;
      }
      content += key_tree;
      $(container).append(content);
      // to display the first object in the file, uncomment the following line
      // setTimeout( function() { showObject(keys[0]['name'],keys[0]['cycle']); }, 20 );
   };

   JSROOTPainter.addDirectoryKeys = function(keys, container, dir_id) {
      var pattern_th1 = /TH1/g;
      var pattern_th2 = /TH2/g;
      var tree_link = '';
      var content = "<p><a href='javascript: key_tree.openAll();'>open all</a> | <a href='javascript: key_tree.closeAll();'>close all</a></p>";
      var k = key_tree.aNodes.length;
      var dir_name = key_tree.aNodes[dir_id]['title'];
      for (var i=0; i<keys.length; ++i) {
         var disp_name = keys[i]['name'];
         keys[i]['name'] = dir_name + '/' + keys[i]['name'];
         var message = keys[i]['className']+' is not yet implemented.';
         tree_link = "javascript:  alert('" + message + "')";
         var node_img = source_dir+'img/page.gif';
         var node_title = keys[i]['className'];
         if (keys[i]['className'].match(/\bTH1/) ||
             keys[i]['className'].match(/\bRooHist/)) {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/histo.png';
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'].match(/\bTH2/)) {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/histo2d.png';
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'].match(/\bTH3/)) {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/histo3d.png';
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'].match(/\bTGraph/) ||
             keys[i]['className'].match(/\RooCurve/)) {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/graph.png';
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'] ==  'TProfile') {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/profile.png';
            node_title = keys[i]['name'];
         }
         else if (keys[i]['name'] == 'StreamerInfo') {
            tree_link = "javascript: displayStreamerInfos(gFile.fStreamerInfo.fStreamerInfos);";
            node_img = source_dir+'img/question.gif';
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'] == 'TDirectory' || keys[i]['className'] == 'TDirectoryFile') {
            tree_link = "javascript: showDirectory('"+keys[i]['name']+"',"+keys[i]['cycle']+","+k+");";
            node_img = source_dir+'img/folder.gif';
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'] == 'TList' || keys[i]['className'] == 'TObjArray') {
            tree_link = "javascript: showCollection('"+keys[i]['name']+"',"+keys[i]['cycle']+","+k+");";
            node_img = source_dir+'img/folder.gif';
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'] == 'TTree' || keys[i]['className'] == 'TNtuple') {
            tree_link = "javascript: readTree('"+keys[i]['name']+"',"+keys[i]['cycle']+","+k+");";
            node_img = source_dir+'img/tree.png';
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'].match('TCanvas')) {
            tree_link = "javascript: showObject('"+keys[i]['name']+"',"+keys[i]['cycle']+");";
            node_img = source_dir+'img/canvas.png';
            node_title = keys[i]['name'];
         }
         if (keys[i]['name'] != '' && keys[i]['className'] != 'TFile') {
            if (keys[i]['className'] == 'TDirectory' || keys[i]['className'] == 'TList' ||
                keys[i]['className'] == 'TObjArray' || keys[i]['className'] == 'TDirectoryFile')
               key_tree.add(k, dir_id, disp_name+';'+keys[i]['cycle'], tree_link, node_title, '', node_img,
                            source_dir+'img/folderopen.gif');
            else if (keys[i]['className'] == 'TNtuple' || keys[i]['className'] == 'TTree')
               key_tree.add(k, dir_id, disp_name+';'+keys[i]['cycle'], tree_link, node_title, '', node_img, node_img);
            else
               key_tree.add(k, dir_id, disp_name+';'+keys[i]['cycle'], tree_link, node_title, '', node_img);
            k++;
         }
      }
      content += key_tree;
      $(container).append(content);
      key_tree.openTo(dir_id, true);
   };

   JSROOTPainter.addCollectionContents = function(list, container, dir_id) {
      var pattern_th1 = /TH1/g;
      var pattern_th2 = /TH2/g;
      var tree_link = '';
      var content = "<p><a href='javascript: key_tree.openAll();'>open all</a> | <a href='javascript: key_tree.closeAll();'>close all</a></p>";
      var k = key_tree.aNodes.length;
      var dir_name = key_tree.aNodes[dir_id]['title'];
      for (var i=0; i<list.length; ++i) {
         var disp_name = list[i]['fName'];
         list[i]['_name'] = dir_name + '/' + list[i]['fName'];
         var message = list[i]['_typename']+' is not yet implemented.';
         tree_link = "javascript:  alert('" + message + "')";
         var node_img = source_dir+'img/page.gif';
         var node_title = list[i]['_typename'];
         if (list[i]['_typename'].match(/\bTH1/) ||
             list[i]['_typename'].match(/\bTH2/) ||
             list[i]['_typename'].match(/\bTH3/) ||
             list[i]['_typename'].match(/\bTGraph/) ||
             list[i]['_typename'].match(/\bRooHist/) ||
             list[i]['_typename'].match(/\RooCurve/)) {
            //tree_link = "javascript: displayMappedObject('"+list[i]['_name']+"');";
            tree_link = "javascript: displayMappedObject('"+list[i]['_name']+"','"+list[i]['_listname']+"',"+list[i]['pos']+");";
            node_img = source_dir+'img/graph.png';
            node_title = list[i]['_name'];
         }
         else if (list[i]['_typename'] == 'TList' || list[i]['_typename'] == 'TObjArray') {
            tree_link = "javascript: showCollection('"+list[i]['_name']+"', 0, "+k+");";
            node_img = source_dir+'img/folder.gif';
            node_title = list[i]['_name'];
         }
         if (list[i]['fName'] != '' && list[i]['_typename'] != 'TFile') {
            key_tree.add(k, dir_id, disp_name, tree_link, node_title, '', node_img);
            k++;
         }
      }
      content += key_tree;
      $(container).append(content);
      key_tree.openTo(dir_id, true);
   };

   JSROOTPainter.displayStreamerInfos = function(streamerInfo, container) {

      delete d;
      var content = "<p><a href='javascript: d_tree.openAll();'>open all</a> | <a href='javascript: d_tree.closeAll();'>close all</a></p>";
      d_tree = new dTree('d_tree');
      d_tree.config.useCookies = false;
      d_tree.add(0, -1, 'Streamer Infos');

      var k = 1;
      var pid = 0;
      var cid = 0;
      var key;
      for (key in streamerInfo) {
         var entry = streamerInfo[key]['name'];
         d_tree.add(k, 0, entry); k++;
      }
      var j=0;
      for (key in streamerInfo) {
         if (typeof(streamerInfo[key]['checksum']) != 'undefined')
            d_tree.add(k, j+1, 'Checksum: ' + streamerInfo[key]['checksum']); ++k;
         if (typeof(streamerInfo[key]['classversion']) != 'undefined')
            d_tree.add(k, j+1, 'Class Version: ' + streamerInfo[key]['classversion']); ++k;
         if (typeof(streamerInfo[key]['title']) != 'undefined')
            d_tree.add(k, j+1, 'Title: ' + streamerInfo[key]['title']); ++k;
         if (typeof(streamerInfo[key]['elements']) != 'undefined') {
            d_tree.add(k, j+1, 'Elements'); pid=k; ++k;
            for (var l=0; l<streamerInfo[key]['elements']['array'].length; ++l) {
               if (typeof(streamerInfo[key]['elements']['array'][l]['element']) != 'undefined') {
                  d_tree.add(k, pid, streamerInfo[key]['elements']['array'][l]['element']['name']); cid=k; ++k;
                  d_tree.add(k, cid, streamerInfo[key]['elements']['array'][l]['element']['title']); ++k;
                  d_tree.add(k, cid, streamerInfo[key]['elements']['array'][l]['element']['typename']); ++k;
               }
               else if (typeof(streamerInfo[key]['elements']['array'][l]['name']) != 'undefined') {
                  d_tree.add(k, pid, streamerInfo[key]['elements']['array'][l]['name']); cid=k; ++k;
                  d_tree.add(k, cid, streamerInfo[key]['elements']['array'][l]['title']); ++k;
                  d_tree.add(k, cid, streamerInfo[key]['elements']['array'][l]['typename']); ++k;
               }
            }
         }
         else if (typeof(streamerInfo[key]['array']) != 'undefined') {
            for (var l=0; l<streamerInfo[key]['array'].length; ++l) {
               d_tree.add(k, j+1, streamerInfo[key]['array'][l]['str']); ++k;
            }
         }
         ++j;
      }
      content += d_tree;
      $(container).html(content);
   };

   var style = "<style>\n"
      +".xaxis path, .xaxis line, .yaxis path, .yaxis line, .zaxis path, .zaxis line {\n"
      +"   fill: none;\n"
      +"   stroke: #000;\n"
      +"   shape-rendering: crispEdges;\n"
      +"}\n"
      +".brush .extent {\n"
      +"  stroke: #fff;\n"
      +"  fill-opacity: .125;\n"
      +"  shape-rendering: crispEdges;\n"
      +"}\n"
      +"rect.zoom {\n"
      +"  stroke: steelblue;\n"
      +"  fill-opacity: 0.1;\n"
      +"}\n"
      /* Correct overflow not hidden in IE9 */
      +"svg:not(:root) { overflow: hidden; }\n"
      +"</style>\n";
   $(style).prependTo("body");

})();


// JSROOTD3Painter.js ends


