// Overrides: JSRootD3Painter.js
//No GUI needed
JSROOTPainter.displayListOfKeys = function(keys) {
      delete key_tree;
      key_tree = new dTree('key_tree');
      key_tree.config.useCookies = false;
      key_tree.add(0, -1, 'File Content');
      var k = 1;
      var tree_link = '';
      for (var i=0; i<keys.length; ++i) {
         if (keys[i]['className'] ==  'TF1') {
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'] ==  'TProfile') {
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'].match('TCanvas')) {
            node_title = keys[i]['name'];
         }
         if (keys[i]['name'] != '' && keys[i]['className'] != 'TFile')
            key_tree.add(k, 0, keys[i]['name']+';'+keys[i]['cycle'], tree_link, keys[i]['name'], '');
         k++;
      }
};

//No GUI needed
JSROOTPainter.addDirectoryKeys = function(keys, dir_id) {
      var pattern_th1 = /TH1/g;
      var pattern_th2 = /TH2/g;
      var tree_link = '';
      var k = key_tree.aNodes.length;
      var dir_name = key_tree.aNodes[dir_id]['title'];
      for (var i=0; i<keys.length; ++i) {
         var disp_name = keys[i]['name'];
         keys[i]['name'] = dir_name + '/' + keys[i]['name'];
         var node_title = keys[i]['className'];
         if (keys[i]['className'].match(/\bTH1/) ||
             keys[i]['className'].match(/\bRooHist/)) {
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'].match(/\bTH2/)) {
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'].match(/\bTH3/)) {
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'].match(/\bTGraph/) ||
             keys[i]['className'].match(/\RooCurve/)) {
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'] ==  'TProfile') {
            node_title = keys[i]['name'];
         }
         else if (keys[i]['name'] == 'StreamerInfo') {
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'] == 'TDirectory' || keys[i]['className'] == 'TDirectoryFile') {
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'] == 'TList' || keys[i]['className'] == 'TObjArray') {
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'] == 'TTree' || keys[i]['className'] == 'TNtuple') {
            node_title = keys[i]['name'];
         }
         else if (keys[i]['className'].match('TCanvas')) {
            node_title = keys[i]['name'];
         }
         if (keys[i]['name'] != '' && keys[i]['className'] != 'TFile') {
            key_tree.add(k, dir_id, disp_name+';'+keys[i]['cycle'], tree_link, node_title, '');
         k++;
         }
      }
      readHistograms();
   };

JSROOTPainter.drawThreeObject = function(obj, idx) {
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
         svg = JSROOTPainter.createCanvas($(render_to), idx);

         if (svg == null) return false;
         JSROOTPainter.drawHStack(svg, null, obj, null);
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

// Overrides: dtree.js
//No GUI needed
dTree.prototype.add = function(id, pid, name, url, title, target) {
   this.aNodes[this.aNodes.length] = new Node(id, pid, name, url, title, target);
};

//Overrides: JSRootIOEvolution.js
//Creating a stack object to draw
function initStackObject(obj_name){
   obj = {};
   obj['fHists'] = [];
   obj['fMaximum'] = -1111;
   obj['fMinimum'] = -1111;
   obj['fTitle'] = obj_name;
   obj['_typename'] = "JSROOTIO.THStack";//JSROOTIO.THStack & JSROOTIO.TH1
   JSROOTCore.addMethods(obj);
   obj['buildStack'] = function() {
      //  build sum of all histograms
      //  Build a separate list fStack containing the running sum of all histograms
      if ('fStack' in this) return;
      if (!'fHists' in this) return;
      var nhists = this['fHists'].length;
      if (nhists <= 0) return;
      this['fStack'] = new Array();
      var h = JSROOTCore.clone(this['fHists'][0]);
      this['fStack'].push(h);
      for (var i=1;i<nhists;i++) {
         h = JSROOTCore.clone(this['fHists'][i]);
         h.add(this['fStack'][i-1]);
         this['fStack'].splice(i, 1, h);
      }
   }
   return obj;
}
//                Red/   Blue  /Green
//               Down/ Nominal /Up
var stackColors = [632, 600, 416];

JSROOTIO.RootFile.prototype.ReadThreeObject = function(obj_name, nuissance, cycle, node_id) {
      // read any object from a root file
      var key = this.GetKey(obj_name+"Down", cycle);
      var key2 = this.GetKey(obj_name.replace(nuissance,''), cycle);
      var key3 = this.GetKey(obj_name+"Up", cycle);
      var colorIndex = 0;
      var stack = {};
      stack = initStackObject(obj_name);
      if (key == null || key2 == null || key3 == null) return;
      this.fTagOffset = key.keyLen;
      var callback = function(file, objbuf) {
         if (objbuf && objbuf['unzipdata']) {
            if (key['className'] == 'TCanvas') {
               var canvas = JSROOTIO.ReadTCanvas(objbuf['unzipdata'], 0);
               if (canvas && canvas['fPrimitives']) {
                  if(canvas['fName'] == "")
                     canvas['fName'] = obj_name;
                  displayThreeObject(canvas, cycle, obj_index);
                  obj_list.push(obj_name+cycle);
                  obj_index++;
               }
            }
            else if (JSROOTIO.GetStreamer(key['className'])) {
               var obj = {};
               obj['_typename'] = 'JSROOTIO.' + key['className'];
               JSROOTIO.GetStreamer(key['className']).Stream(obj, objbuf['unzipdata'], 0);
               if (key['className'] == 'TFormula') {
                  JSROOTCore.addFormula(obj);
               }
               else if (key['className'] == 'TNtuple' || key['className'] == 'TTree') {
                  displayTree(obj, cycle, node_id);
               }
               else {
                  JSROOTCore.addMethods(obj);
                  obj['fLineColor'] = stackColors[colorIndex];
                  stack['fHists'].push(obj);
                  if (stack['fHists'].length == 3){
                     displayThreeObject(stack, cycle, obj_index);
                     obj_list.push(obj_name+cycle);
                     obj_index++;
                  }
                  colorIndex++;
               }
            }
            delete objbuf['unzipdata'];
            objbuf['unzipdata'] = null;
         }
      };
      this.ReadObjBuffer(key, callback);
      this.fTagOffset = key2.keyLen;
      this.ReadObjBuffer(key2, callback);
      this.fTagOffset = key3.keyLen;
      this.ReadObjBuffer(key3, callback);
};
//Multiple titles in multiple histograms, need only 1
stackTitleCounter = 0;
JSROOTPainter.drawTitle = function(vis, histo, pad) {
   /* draw the title only if we don't draw from a pad (see Olivier for details) */
   var w = vis.attr("width"), h = vis.attr("height");
   var font_size = Math.round(0.050 * h);
   var l_title = this.translateLaTeX(histo['fTitle']);

   if (stackTitleCounter == 0){
      if (!pad || typeof(pad) == 'undefined') {
         vis.append("text")
            .attr("class", "title")
            .attr("text-anchor", "middle")
            .attr("x", w/2)
            .attr("y", 1 + font_size)
            .attr("font-family", "Arial")
            .attr("font-size", font_size)
            .text(l_title.replace("Down", ''));
      }
      stackTitleCounter++;
   }else if(stackTitleCounter == 3){
      stackTitleCounter = 0;
   }
   else
      stackTitleCounter++;
}

//Always draw the same way-> option.Error = false; 
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
      option.Error = false;
      return option;
   };

JSROOTPainter.createCanvas = function(element, idx) {
      var w = element.width(), h = w * 0.6666666;
      var render_to = '#histogram' + idx;
      d3.select(render_to).style("background-color", 'white');
      d3.select(render_to).style("width", "100%");

      var svg = d3.select(render_to).append("svg")
                  .attr("width", w-60)
                  .attr("height", h)
                  .style("background-color", 'white');
      defs = svg.append('svg:defs');
      return svg;
   };