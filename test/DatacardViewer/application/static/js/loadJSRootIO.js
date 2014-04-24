var gFile;
var obj_list = new Array();
var obj_index = 0;
var last_index = 0;
var function_list = new Array();
var func_list = new Array();
var frame_id = 0;
var random_id = 0;
var source_dir = "static/lib/JSRootIO-2.1/";

function loadScript(url, callback){
   // dynamic script loader using callback
   // (as loading scripts may be asynchronous)
   var script = document.createElement("script");
   script.type = "text/javascript";
   if (script.readyState){ // Internet Explorer specific
      script.onreadystatechange = function(){
         if (script.readyState == "loaded" || script.readyState == "complete"){
            script.onreadystatechange = null;
            callback();
         }
      };
   }else{ // Other browsers
      script.onload = function(){
         callback();
      };
   }
   var rnd = Math.floor(Math.random()*80000);
   script.src = url;
   document.getElementsByTagName("head")[0].appendChild(script);
}

function assertPrerequisitesAndRead(){
   if (typeof JSROOTIO == "undefined"){
      // if JSROOTIO is not defined, then dynamically load the required scripts
      loadScript('http://d3js.org/d3.v2.min.js', function(){
         loadScript(source_dir+'jquery.mousewheel.js', function(){
            loadScript(source_dir+'dtree.js', function(){
               loadScript(source_dir+'rawinflate.js', function(){
                  loadScript(source_dir+'JSRootCore.js', function(){
                     loadScript(source_dir+'three.min.js', function(){
                        loadScript(source_dir+'JSRootD3Painter.js', function(){
                           loadScript(source_dir+'JSRootIOEvolution.js', function(){
                              loadScript('static/js/JSRootIO_Overrides.js', function(){
                                 readRootFiles();
                              });
                           });
                        });
                     });
                  });
               });
            });
         });
      });
   }else{
      readRootFiles();
   }
}

//After Root files keys read callback
function userCallback(file){
   readRootContent(file.fKeys);
}

function readRootFiles(){
   var navigator_version = navigator.appVersion;
   if (typeof ActiveXObject == "function"){ // Windows
      // detect obsolete browsers
      if ((navigator_version.indexOf("MSIE 8") != -1) || 
          (navigator_version.indexOf("MSIE 7") != -1)){
         alert("You need at least MS Internet Explorer version 9.0. Note you can also use any other web browser (except Opera)");
         return;
      }
   }else{
      // Safari 5.1.7 on MacOS X doesn't work properly
      if ((navigator_version.indexOf("Windows NT") == -1) &&
          (navigator_version.indexOf("Safari") != -1) &&
          (navigator_version.indexOf("Version/5.1.7") != -1)){
         alert("There are know issues with Safari 5.1.7 on MacOS X. It may become unresponsive or even hangs. You can use any other web browser (except Opera)");
         return;
      }
   }
   //TODO multiple files
   if (gFile) {
      gFile.Delete();
      delete gFile;
   }
   if (rootjsFiles.length>1)
      alert("Multiple ROOT files no supported");
   else
      for (var i = 0;i<rootjsFiles.length;i++){
         gFile = null;
         gFile = new JSROOTIO.RootFile("datacards/files/"+rootjsFiles[i]);
      }
}

function displayDirectory(directory, cycle, dir_id) {
   JSROOTPainter.addDirectoryKeys(directory.fKeys, dir_id);
}

function showDirectory(dir_name, cycle, dir_id) {
   gFile.ReadDirectory(dir_name, cycle, dir_id);
}

function readRootContent(keys){
   if(keys[0]["className"] == "TDirectoryFile" || keys[0]["className"] == "TDirectory")
      for (aBin in datacardShapeMap){
         for (var i = 0;i<keys.length;i++){
            if (aBin == keys[i]["name"]){
               showDirectory(aBin, 1, i+1);
            }
         }
      }
   else
      readHistograms();
}

function getHistogramNumber(binProcNuisArr){
   var index = 0;
   var histoBin = binProcNuisArr[0];
   var histoProc = binProcNuisArr[1];
   var histoNuis = binProcNuisArr[2];
   for (aBin in datacardShapeMap){
      for (aProcess in datacardShapeMap[aBin]){
         var tempNuis = datacardShapeMap[aBin][aProcess].slice(3);
         for (var j = 0; j<tempNuis.length; j++){
            if (aBin == histoBin && aProcess == histoProc && tempNuis[j] == histoNuis){
               return index;
            }
            else{
               index++;
            }
         }
      }
   }
}

function getHistogramPath(aBin, aProcess, aNuissance){
   var path = datacardShapeMap[aBin][aProcess][2];
   path = path.replace("$CHANNEL", aBin);
   path = path.replace("$PROCESS", aProcess);
   path = path.replace("$MASS", settings_mass);
   path = path.replace("$SYSTEMATIC", aNuissance);
   return path;
}
//TODO parse better datacardShapeMap, add Up/Down/Nominal combined
function readHistograms(){
   for (aBin in datacardShapeMap){
      for (aProcess in datacardShapeMap[aBin]){
         var tempNuis = datacardShapeMap[aBin][aProcess].slice(3);
         for (var j = 0; j<tempNuis.length; j++){
            showObject(getHistogramPath(aBin, aProcess, tempNuis[j])+"Up",1);
         }
      }
   }
}

function displayListOfKeys(keys) {
   JSROOTPainter.displayListOfKeys(keys);
}

function findObject(obj_name) {
   for (var i in obj_list) {
      if (obj_list[i] == obj_name) {
         var findElement = $('#report').find('#histogram'+i);
         if (findElement.length) {
            var element = findElement[0].previousElementSibling.id;
            showElement('#'+element);
            return true;
         }
      }
   }
   return false;
}

function showElement(element) {
   if ($(element).next().is(":hidden")) {
      $(element)
         .toggleClass("ui-accordion-header-active ui-state-active ui-state-default ui-corner-bottom")
         .find("> .ui-icon").toggleClass("ui-icon-triangle-1-e ui-icon-triangle-1-s").end()
         .next().toggleClass("ui-accordion-content-active").slideDown(0);
   }
   $(element)[0].scrollIntoView();
}

function displayRootStatus(msg) {
   $("#status").append(msg);
};

function displayStreamerInfos(streamerInfo) {
   var findElement = $('#report').find('#treeview');
   if (findElement.length) {
      var element = findElement[0].parentElement.previousSibling.id;
      showElement('#'+element);
   }
   else {
      var uid = "uid_accordion_"+(++last_index);
      var entryInfo = "<h5 id=\""+uid+"\"><a> Streamer Infos </a>&nbsp; </h5><div>\n";
      entryInfo += "<h6>Streamer Infos</h6><span id='treeview' class='dtree'></span></div>\n";
      $("#report").append(entryInfo);
      JSROOTPainter.displayStreamerInfos(streamerInfo, '#treeview');
   }
};

function findObject(obj_name) {
   for (var i in obj_list) {
      if (obj_list[i] == obj_name) {
         var findElement = $('#report').find('#histogram'+i);
         if (findElement.length) {
            var element = findElement[0].previousElementSibling.id;
            showElement('#'+element);
            return true;
         }
      }
   }
   return false;
};

function showObject(obj_name, cycle) {
   gFile.ReadObject(obj_name, cycle);
};

function displayCollection(cont, cycle, c_id) {
   var url = $("#urlToLoad").val();
   $("#status").html("file: " + url + "<br/>");
   JSROOTPainter.addCollectionContents(cont, '#status', c_id);
};

function showCollection(name, cycle, id) {
   gFile.ReadCollection(name, cycle, id);
};

function readTree(tree_name, cycle, node_id) {
   gFile.ReadObject(tree_name, cycle, node_id);
};

function displayTree(tree, cycle, node_id) {
   var url = $("#urlToLoad").val();
   $("#status").html("file: " + url + "<br/>");
   JSROOTPainter.displayTree(tree, '#status', node_id);
};

function displayObject(obj, cycle, idx) {
   if (!obj['_typename'].match(/\bJSROOTIO.TH1/) &&
       !obj['_typename'].match(/\bJSROOTIO.TH2/) &&
       !obj['_typename'].match(/\bJSROOTIO.TH3/) &&
       !obj['_typename'].match(/\bJSROOTIO.TGraph/) &&
       !obj['_typename'].match(/\bRooHist/) &&
       !obj['_typename'].match(/\RooCurve/) &&
       obj['_typename'] != 'JSROOTIO.TCanvas' &&
       obj['_typename'] != 'JSROOTIO.TF1' &&
       obj['_typename'] != 'JSROOTIO.TProfile') {
      if (typeof(checkUserTypes) != 'function' || checkUserTypes(obj) == false)
         return;
   }
   var entryInfo = "<div id='histogram" + idx + "'>\n";
   $("#report").append(entryInfo);
   JSROOTPainter.drawObject(obj, idx);
   $("#histogram" + idx).hide();
};

function displayMappedObject(obj_name, list_name, offset) {
   var obj = null;
   for (var i=0; i<gFile['fObjectMap'].length; ++i) {
      if (gFile['fObjectMap'][i]['obj']['_name'] == obj_name) {
         obj = gFile['fObjectMap'][i]['obj'];
         break;
      }
   }
   if (obj == null) {
      gFile.ReadCollectionElement(list_name, obj_name, 1, offset);
      return;
   }
   if (!obj['_typename'].match(/\bJSROOTIO.TH1/) &&
       !obj['_typename'].match(/\bJSROOTIO.TH2/) &&
       !obj['_typename'].match(/\bJSROOTIO.TH3/) &&
       !obj['_typename'].match(/\bJSROOTIO.TGraph/) &&
       !obj['_typename'].match(/\bRooHist/) &&
       !obj['_typename'].match(/\RooCurve/) &&
       obj['_typename'] != 'JSROOTIO.TCanvas' &&
       obj['_typename'] != 'JSROOTIO.TF1' &&
       obj['_typename'] != 'JSROOTIO.TProfile') {
      if (typeof(checkUserTypes) != 'function' || checkUserTypes(obj) == false)
         return;
   }
   var uid = "uid_accordion_"+(++last_index);
   var entryInfo = "<h5 id=\""+uid+"\"><a> " + obj['fName'] + "</a>&nbsp; </h5>\n";
   entryInfo += "<div id='histogram" + obj_index + "'>\n";
   $("#report").append(entryInfo);
   //draw
   JSROOTPainter.drawObject(obj, obj_index);
   obj_list.push(obj['fName']);
   obj_index++;
};