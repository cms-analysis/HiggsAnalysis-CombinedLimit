// Overrides: JSRootD3Painter.js
//
// core methods overriten for Javascript ROOT Graphics, using d3.js.
//

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
//

dTree.prototype.add = function(id, pid, name, url, title, target) {
   this.aNodes[this.aNodes.length] = new Node(id, pid, name, url, title, target);
};