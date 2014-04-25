// Overrides: JSRootD3Painter.js

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

dTree.prototype.add = function(id, pid, name, url, title, target) {
   this.aNodes[this.aNodes.length] = new Node(id, pid, name, url, title, target);
};

//Overrides: JSRootIOEvolution.js
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
//              Red/     Blue /Green
//              Nominal/ Up   /Down
var stackColors = [632, 600, 416];

JSROOTIO.RootFile.prototype.ReadThreeObject = function(obj_name, nuissance, cycle, node_id) {
      // read any object from a root file
      var key = this.GetKey(obj_name.replace(nuissance,''), cycle);
      var key2 = this.GetKey(obj_name+"Up", cycle);
      var key3 = this.GetKey(obj_name+"Down", cycle);
      var colorIndex = 0;
      var stack = {};
      stack = initStackObject(obj_name);
      //stack['fOption'] =
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

//options.Error = false; to fix in JSROOTIOEVO
stackTitleCounter = 0;
JSROOTPainter.drawTitle = function(vis, histo, pad) {
   /* draw the title only if we don't draw from a pad (see Olivier for details) */
   var w = vis.attr("width"), h = vis.attr("height");
   var font_size = Math.round(0.050 * h);
   var l_title = this.translateLaTeX(histo['fTitle']);

   if (stackTitleCounter == 1){
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
