
var $datacardTable;
var bin;
var observation;
var bins;
var processes;
var signals;
var rate;
var nuisances;
var shapeMap;
var colors;
var colorsToShow;
var colIDforDefault;
var colIDtoShow;
var colValuestoShow;
var colHeaderstoShow;
var rowHeaders;
var rowHeadersMaxWidth;
var datacardShapeMap;
var rootjsFiles;
var settings_mass = 127;//120.7;
var histograms;
var histogramsWidth = {};

/*
   Function: init_variables
   Initializes global variables for one datacard.
*/
function init_variables(){
    $datacardTable = $('#datacard');
    //Parced main datacard variables
    bin = [];
    observation = [];
    bins = [];
    processes = [];
    signals = [];
    rate = [];
    nuisances = [];
    shapeMap = [];
    //Datacard variables for table functionality 
    colIDforDefault = [];
    colIDtoShow = [];
    colValuestoShow = [];
    colHeaderstoShow =[];
    rowHeaders = [];
    rowHeadersMaxWidth = 0;
    datacardShapeMap = {};
    rootjsFiles = [];
    histograms = {};
    //Histogram plots numbering
    $("#report").empty();
    if(typeof obj_index !== 'undefined'){
        obj_index = 0;
    };
    
}

/*
   Function: get_colors

   Associates datacard's processes signals with colors. 

   Parameters:

      signals - Array of signal integers.

   Returns:

      An array of colors for processes.

*/
function get_colors(signals){
    var backgrounds = ["#DDA0DD", "#FCF8E3", "#DFF0D8", "#D9EDF7", "#F2DEDE", "#FFE4C4",
                        "#DEB887", "#BDB76B", "#E9967A", "#8FBC8F", "#CD5C5C", "F08080",
                        "#20B2AA", "#778899", "#191970", "#FF4500", "#D8BFD8", "#708090",
                        'rgb(122, 142, 153)', 'rgb(104, 130, 150)', 'rgb(109, 122, 132)',
                        'rgb(124, 153, 209)', 'rgb(127, 127, 155)', 'rgb(170, 165, 191)',
                        'rgb(211, 206, 135)', 'rgb(221, 186, 135)', 'rgb(188, 158, 130)',
                        'rgb(198, 153, 124)','rgb(191, 130, 119)', 'rgb(206, 94, 96)',
                        'rgb(170, 142, 147)','rgb(165, 119, 122)', 'rgb(147, 104, 112)'
    ];
    var colors = [];
    for(var i=0;i<signals.length;i++)
        colors.push(backgrounds[signals[i]]);
    return colors;
}

/*
   Function: eliminateDuplicates

   Eliminates duplicate values from array. 

   Parameters:

      arr - Array of values.

   Returns:

      An array of values without duplicates.

*/
function eliminateDuplicates(arr){
    obj ={};
    out =[];
    for (var i=0;i<arr.length;i++) {
        obj[arr[i]]=0;
    }
    for (var i in obj) {
        out.push(i);
    }
    return out;
}

/*
   Function: countDuplicates

   Counts duplicate values in array. 

   Parameters:

      arr - Array of values.

   Returns:

      An object of [values] = count.

*/
function countDuplicates(arr){
    obj ={};
    out =[];
    for (var i=0;i<arr.length;i++) {
        if (typeof obj[arr[i]] === 'undefined') 
            obj[arr[i]]=1;
        else
            obj[arr[i]]=obj[arr[i]]+1;
    }
    for (var i=0;i<bin.length;i++)
        if (typeof obj[bin[i]] !== 'undefined') 
            out.push(obj[bin[i]]);
    return out;
}

/*
   Function: getColDataIndicesOf

   Searches array for a given value and saves that value index in array. 
   This function is used for dynamic handsontable columns changing.
   
   Parameters:

      searchArr - value.
	  arr - Array of values.

   Returns:

      An array of objects ["data"] = position.
	  

*/
function getColDataIndicesOf(searchArr, arr) {
    var indices = [];
    for (var i = 0; i< arr.length; i++){
        if (arr[i]==searchArr)
            indices.push({"data":i});
    }
    return indices;
}

/*
   Function: getColDataIndicesOf

   Searches array for a given value and saves that value index in array.  

   Parameters:

      searchArr - value.
	  arr - Array of values.

   Returns:

      An array of positions.

*/
function getColIndicesOf(searchArr, arr) {
    var indices = [];
    for (var i = 0; i< arr.length; i++){
        if (arr[i]==searchArr)
            indices.push(i);
    }
    return indices;
}

$.fn.textWidth = function(){
    var html_org = $(this).html();
    var html_calc = '<span>' + html_org + '</span>';
    $(this).html(html_calc);
    var width = $(this).find('span:first').width();
    $(this).html(html_org);
    return width;
};

String.prototype.insert = function (index, string) {
  if (index > 0)
    return this.substring(0, index) + string + this.substring(index, this.length);
  else
    return string + this;
};

/*
   Function: init_JSRootIO
   Initializes JSRootIO from loadJSRootIO.js file.  
*/
function init_JSRootIO(){
    assertPrerequisitesAndRead();
    //TODO mass settings
    //set_settings();
}

/*
   Function: show_datacard

   Main function which initializes, parses and generates datacard.  

   Parameters:

      data - a datacard object returned from the server.
*/
function show_datacard(data){
    init_variables();
    parse_datacard(data);
    generate_datacard();
}

/*
   Function: parse_datacard

   Parses the data from the server into global variables.  

   Parameters:

      data - a datacard object returned from the server.
*/
function parse_datacard(data){
    $("h2#filename").html(data.filename);
    for (var abin in data.binsProcessesRates){
        bin.push(abin);
        observation.push(data.observation[abin]);
        var total = 0;
        for(var i in data.binsProcessesRates[abin]) {
            if (data.binsProcessesRates[abin].hasOwnProperty(i)) {
                total++; 
            }
        }
        var signalNr = 1;
        for (var aProcess in data.binsProcessesRates[abin]){ 
            bins.push(abin);
            processes.push(aProcess);
            rate.push(data.binsProcessesRates[abin][aProcess]);
            if (signalNr>total-1){
                signals.push(0);
            }
            else{
                signals.push(signalNr);
                signalNr++;
            }
        }
    }
    var shapesIDs = [];
    for (var i = 0; i<data.nuisances.length; i++){
        var temp = [];
        temp.push(data.nuisances[i][0]);
        temp.push(data.nuisances[i][2]);
        var index = 0;
        for (var aBin in data.nuisances[i][4]){
            for (var aProcess in data.nuisances[i][4][aBin]){
                if (data.nuisances[i][4][aBin][aProcess] == 0)
                    temp.push("-");
                else if(data.nuisances[i][2] === "shape"|| data.nuisances[i][2] === "shapeN2"){
                    temp.push("<button class='btn btn-default btn-xs fa fa-eye' id='"+aBin+":"+aProcess+":"+data.nuisances[i][0]+"'></button>");
                    if (shapesIDs.indexOf(index)<0)
                        shapesIDs.push(index);
                }
                else if((data.nuisances[i][4][aBin][aProcess]) instanceof Array){
                    temp.push(data.nuisances[i][4][aBin][aProcess][0]+"/"+data.nuisances[i][4][aBin][aProcess][1]);
                }
                else
                    temp.push(data.nuisances[i][4][aBin][aProcess]);
                index++;
            }
        }
        nuisances.push(temp);
    }
    //only get the needed shapes to draw
    var shapesToDraw = [];
    for(var i = 0;i<shapesIDs.length;i++)
        if (shapesToDraw[bins[shapesIDs[i]]] === undefined){
            var temp = [];
            temp.push(processes[shapesIDs[i]]);
            shapesToDraw[bins[shapesIDs[i]]] = temp; 
        }else
            shapesToDraw[bins[shapesIDs[i]]].push(processes[shapesIDs[i]]);
    //shapes process/bin/file/histogram-name/histogram-name-for-systematics
    if(Object.getOwnPropertyNames(data.shapeMap).length > 0){
        var shapeBin = sort_obj_keys(data.shapeMap);
        for (var j = 0;j<shapeBin.length;j++){
            //* - rule applies to all processes, unless a more specific rule exists for it
            if (shapeBin[j] == "*")
                for (shapeBinDraw in shapesToDraw)
                    datacardShapeMap[shapeBinDraw]={};
            else
                datacardShapeMap[shapeBin[j]] = {};
            var shapeProc = sort_obj_keys(data.shapeMap[shapeBin[j]]);
            for (var k = 0;k<shapeProc.length;k++){
                //* - rule applies to all channels, unless a more specific rule exists for it
                if (shapeBin[j] == "*" && shapeProc[k] == "*")
                    for (shapeBinDraw in shapesToDraw)
                        for (var i = 0;i<shapesToDraw[shapeBinDraw].length; i++)
                            datacardShapeMap[shapeBinDraw][shapesToDraw[shapeBinDraw][i]] = data.shapeMap[shapeBin[j]][shapeProc[k]].slice(0);
                else if (shapeProc[k] == "*")
                    for (var i = 0;i<shapesToDraw[shapeBin[j]].length; i++)
                        datacardShapeMap[shapeBin[j]][shapesToDraw[shapeBin[j]][i]] = data.shapeMap[shapeBin[j]][shapeProc[k]].slice(0);
                else if(shapeBin[j] == "*") 
                    for (shapeBinDraw in shapesToDraw)
                        datacardShapeMap[shapeBinDraw][shapeProc[k]] = data.shapeMap[shapeBin[j]][shapeProc[k]].slice(0);
                else if (shapeProc[k] == "data_obs")
                    continue;
                else
                    datacardShapeMap[shapeBin[j]][shapeProc[k]] = data.shapeMap[shapeBin[j]][shapeProc[k]].slice(0);

                if (rootjsFiles.indexOf(data.shapeMap[shapeBin[j]][shapeProc[k]][0]) < 0)
                    rootjsFiles.push(data.shapeMap[shapeBin[j]][shapeProc[k]][0]);
            }
        }
        //push nuisances to datacardShapeMap
        for (var i = 0; i<data.nuisances.length; i++){
            if(data.nuisances[i][2] === "shape" || data.nuisances[i][2] === "shapeN2"){
                for (aBin in data.nuisances[i][4]){
                    for (aProcess in data.nuisances[i][4][aBin]){
                        if (data.nuisances[i][4][aBin][aProcess] != 0){
                            datacardShapeMap[aBin][aProcess].push(data.nuisances[i][0]);
                        }
                    }
                }
            }
        }
        init_JSRootIO();
    }
}

/*
   Function: sort_obj_keys

   Sorts object keys and pushes them into array. 

   Parameters:

      obj - Object.

   Returns:

      An array of sorted object keys.

*/
function sort_obj_keys(obj){
    var keys = [];
    for (k in obj){
        if (obj.hasOwnProperty(k)){
            keys.push(k);
        }
    }
    keys.sort();
    return keys;
}

/*
   Function: generate_datacard
   Creates a handsontable with needed functionality (searching, grouping, sorting and menus)
*/
function generate_datacard(){
    var tableData = [];
    var tableWidths = [];
    var columns = [];
    var groups = [];
    colors = get_colors(signals);
    colorsToShow = colors;
    for (var i = 0; i<bins.length; i++){
        tableWidths.push(80);
        columns.push({readOnly: true});
        colIDforDefault.push({"data": i});
    }
    rowHeaders.push("process");
    tableData.push(processes);
    var temp;
    var temp2;
    var temp3;
    for (var i = 0; i<nuisances.length; i++){
        if (nuisances[i][1] == "param")
            continue;
        rowHeaders.push("<div class='pull-left'>"+nuisances[i][0]+": "+nuisances[i][1]+'</div><button class="fa fa-bars pull-right"></button>');
        tableData.push(nuisances[i].slice(2));
        temp = [];
        temp2 = [];
        temp3 = [];
        for(var j = 2; j<nuisances[i].length; j++){
            if(nuisances[i][j]!="-"){
                temp.push({"data": j-2});
                temp2.push(nuisances[i][j]);
                temp3.push(bins[j-2]);
            }
        }
        colIDtoShow.push(temp);
        colValuestoShow.push(temp2);
        colHeaderstoShow.push(temp3);
    }
    if (typeof($datacardTable.handsontable('getInstance'))!=='undefined')
        $datacardTable.handsontable('getInstance').destroy();

    $datacardTable.handsontable({
        rowHeaders: rowHeaders,
        colHeaders: bins,
        minCols: bins.length,
        maxCols: bins.length,
        data: tableData,
        colWidths: tableWidths,
        //manualColumnResize: true,
        //columns: columns,
        stretchH: 'all',
        minSpareRows: 0,
        contextMenu: false,
        cells: function (row, col, prop) {
            if (row === 0){
                this.renderer = processRenderer;
            }else
                this.renderer = cellButtonRenderer;
            return {readOnly: true};
        },
        afterRender: function(isForced){
            var div = $("thead tr:first th:first div:first");
            if(div.html() == "&nbsp;")
                div.append("bin");
            if(!isForced){
                var lastOffsetColumn = this.view.wt.lastOffsetColumn;
                change_on_events_proc(true, lastOffsetColumn);
                change_on_events_bin(true, lastOffsetColumn);
            }
        }
        /*,
        currentRowClassName: 'currentRow',
        currentColClassName: 'currentCol',
        autoWrapRow: true*/
    });
    count_headers_width();
    change_headers_width();
    change_on_events_bin(true);
    change_on_events_proc(true);
    change_on_events_nuisances(true);
    add_search_to_table(tableData);
}

/*
   Function: processRenderer
   Dynamically changes the color of all process cells. 
*/
var processRenderer = function (instance, td, row, col, prop, value, cellProperties) {
    Handsontable.renderers.TextRenderer.apply(this, arguments);
    $(td).css({
        background: colorsToShow[col]
    });
};

/*
   Function: processRenderer
   Dynamically changes the color of all process cells. 
*/
var descriptionRenderer = function (instance, td, row, col, prop, value, cellProperties) {
    var escaped = Handsontable.helper.stringify(value);
    escaped = strip_tags(escaped, '<em><b><a>'); 
    td.innerHTML = escaped;
    return td;
};

/*
   Function: cellButtonRenderer
   Adds buttons to shape cells. 
*/
var cellButtonRenderer = function (instance, td, row, col, prop, value, cellProperties) {
    if((value+"").indexOf("<button")>=0){
        td.innerHTML = Handsontable.helper.stringify((value+"").insert(value.indexOf('>'), ' style="width:'+$(td).width()+'px"'));
        add_cell_dialog_event($(td).find('button:first'));
    }else
        td.innerHTML = value;
};

/*
   Function: add_cell_dialog_event
   
   Adds on click event to shape cell button. 
   
   Parameters:

      button - jQuery button object.
   
*/
function add_cell_dialog_event(button){
    button.off();
    button.on("click", function(e){
        build_cell_dialog(button);
        e.stopPropagation();
    });
}

/*
   Function: build_cell_dialog
   
   Shows histogram plot inside Bootstrap Dialog.
   
   Parameters:

      button - jQuery button object.
   
*/
function build_cell_dialog(button){
    var binProc = button.attr('id').split(":");
    var histNr = getHistogramNumber(binProc);
    var $hist = $("#histogram"+histNr).clone(true);
    var width;
    var dialog = new BootstrapDialog({ 
        title: function(){
            if ($hist.children().length === 0 )
                if (histograms[binProc] === undefined)
                    return "Root file not loaded";
                else
                    return "Bin: "+binProc[0]+", proccess: "+binProc[1]+", nuissance: "+binProc[2]+" histograms";
            else
                return "Bin: "+binProc[0]+", proccess: "+binProc[1]+", nuissance: "+binProc[2]+" histograms";
        },
        message: function(dialog){
            if ($hist.children().length === 0 ){
                if (histograms[binProc] === undefined){
                    width = 600;
                    return "Not yet loaded, close it and try again.\nRoot files with \"Combination of\", \"ProcessID\" and \"RooWorkspace\" are not supported";
                }else{
                    width = histogramsWidth[binProc];
                    return histograms[binProc];
                }
            }else{
                if (histograms[binProc] === undefined){
                    histogramsWidth[binProc] = +$("#histogram"+histNr+" svg").attr("width")+40;
                    width = histogramsWidth[binProc];
                    var html = "<b>Down</b>:<span style='color:red;'>red</span>, ";
                    html+="<b>Nominal</b>: <span style='color:blue;'>blue</span>, ";
                    html+="<b>Up</b>: <span style='color:green;'>green</span>.";
                    histograms[binProc] = $("#histogram"+histNr).show().append(html);
                    return histograms[binProc];
                }
                else{
                    width = histogramsWidth[binProc];
                    return histograms[binProc];
                }
            }
        }
    });
    dialog.realize();
    dialog.getModalDialog().css('width', width+'px');
    dialog.open();
}

/*
   Function: count_headers_width
   Counts all nuisances name width in the table.
*/
function count_headers_width(){
    var rowIndex = 0;
    $("th").each(function() {
        if($(this).html().indexOf('<div class="relative">') < 0
        && $(this).html() != "process"){
            count_header_width($(this));
            rowIndex++;
        }
    });
}

/*
   Function: change_headers_width
   Changes datacard table nuisance column width.
*/
function change_headers_width(){
    $(".handsontable col.rowHeader ").css( "width", rowHeadersMaxWidth + 50);
    $datacardTable.handsontable('getInstance').render();
}

/*
   Function: set_processes_colors
   
   Dynamically changes process colors after sorting/grouping.
      
   Parameters:

      cols - Processes colors array.
*/
function set_processes_colors(cols){
    if(typeof(cols)==='undefined')
        colorsToShow = colors;
    else{
       colorsToShow=[];
        for (var i = 0; i<cols.length;i++){
            colorsToShow.push(colors[cols[i].data]);
        }
    } 
}

/*
   Function: change_on_events_bin_proc
   
   Sorts datacard's handsontable by bin -> process.
      
   Parameters:

      firstClick - Boolean variable to identify if handsontable cell was clicked for the first time.
	  binstoShow - Array of bins that have values for selected process.
	  selectedProc - Array of selected process items to show.
*/
function change_on_events_bin_proc(firstClick, binstoShow, selectedProc){
    var colIndex = 0;
    $("th").each(function() {
        if($(this).html().indexOf('<span class="colHeader">') >= 0 
        && jQuery.inArray($(this).children().children().html(), bin) >= 0){
            var binShow = binstoShow;
            var proc = selectedProc;
            var index = colIndex;
            if (firstClick){
                $(this).off();
                $(this).on("click", function(e){
                    set_processes_colors([binShow[index]]);
                    $datacardTable.handsontable('getInstance').updateSettings({              
                        columns: [binShow[index]],
                        colHeaders: proc[index]
                    });
                    change_on_events_bin_proc(false, binShow, proc);
                    change_on_events_nuisances(false);
                    hide_rows([binShow[index]]);
                });
            }
            else{
                $(this).on("click", function(e){
                    show_all_rows();
                    set_processes_colors(binShow);
                    $datacardTable.handsontable('getInstance').updateSettings({              
                        columns: binShow,
                        colHeaders: proc
                    });
                    change_on_events_proc(false);
                    change_on_events_bin_proc(true, binShow, proc);
                    change_on_events_nuisances(false);
                    hide_rows(binShow);    
                });
            }
            colIndex++;
        }
    });
}

/*
   Function: change_on_events_proc_bin
   
   Sorts datacard's handsontable by process -> bin.
      
   Parameters:

      firstClick - Boolean variable to identify if handsontable cell was clicked for the first time.
	  processestoShow - Array of processes that have values for selected bin.
	  selectedBin - Array of selected bin items to show.
*/
function change_on_events_proc_bin(firstClick, processestoShow, selectedBin){
    var procIndex = 0;
    $("table.htCore tbody tr:first td").each(function() {
        var proc = processestoShow;
        var bin = selectedBin;
        var index = procIndex;
        if (firstClick){
            $(this).off();
            $(this).on("click", function(e){
                set_processes_colors([proc[index]]);
                $datacardTable.handsontable('getInstance').updateSettings({              
                    columns: [proc[index]],
                    colHeaders: bin
                });
                change_on_events_proc_bin(false, proc, bin);
                change_on_events_nuisances(false);
                hide_rows([proc[index]]);
            });
        }
        else{
            $(this).off();
            $(this).on("click", function(e){
                show_all_rows();
                set_processes_colors(proc);
                $datacardTable.handsontable('getInstance').updateSettings({              
                    columns: proc,
                    colHeaders: bin
                });
                change_on_events_bin(false);
                change_on_events_proc_bin(true, proc, bin);
                change_on_events_nuisances(false);
                hide_rows(proc);    
            });
        }
        procIndex++;
    });
}

/*
   Function: change_on_events_bin
   
   Sorts datacard's handsontable by bin.
      
   Parameters:

      firstClick - Boolean variable to identify if handsontable cell was clicked for the first time.
	  lastOffsetColumn - Integer of column offset.
*/
function change_on_events_bin(firstClick, lastOffsetColumn){
    var colIndex;
    if (typeof(lastOffsetColumn)!=='undefined')
        colIndex = lastOffsetColumn;
    else
        colIndex = 0;
    $("th").each(function() {
        if($(this).html().indexOf('<span class="colHeader">') >= 0 
        && jQuery.inArray($(this).children().children().html(), bin) >= 0){
            if(!firstClick){
                $(this).off();
                $(this).on("click", function(e){
                    remove_select();                
                    show_all_rows();
                    set_processes_colors();
                    $datacardTable.handsontable('getInstance').updateSettings({               
                        columns: colIDforDefault,
                        rowHeaders: rowHeaders,
                        colHeaders: bins
                    });
                    change_on_events_bin(true);
                    change_on_events_proc(true);
                    change_on_events_nuisances(true);
                });
            }
            else{
                $(this).off();
                var selectedBin = [];
                var processestoShow = [];
                var selected = getColIndicesOf(bins[colIndex], bins);
                for (var i = 0; i<selected.length;i++){
                    selectedBin.push(bins[selected[i]]);
                }
                processestoShow = getColDataIndicesOf(bins[colIndex], bins);
                $(this).on("click", function(e){
                    set_processes_colors(processestoShow);
                    $datacardTable.handsontable('getInstance').updateSettings({              
                        columns: processestoShow,
                        colHeaders: selectedBin
                    });
                    change_on_events_bin(false);
                    change_on_events_proc_bin(true, processestoShow, selectedBin);
                    change_on_events_nuisances(false);
                    hide_rows(processestoShow);    
                });
            }    
            colIndex++;
        }
    });
}

/*
   Function: change_on_events_proc
   
   Sorts datacard's handsontable by process.
      
   Parameters:

      firstClick - Boolean variable to identify if handsontable cell was clicked for the first time.
	  lastOffsetColumn - Integer of column offset.
*/
function change_on_events_proc(firstClick, lastOffsetColumn){
    var procIndex;
    if (typeof(lastOffsetColumn)!=='undefined')
        procIndex = lastOffsetColumn;
    else
        procIndex = 0;
    $("table.htCore tbody tr:first td").each(function() {
        if(!firstClick){
            $(this).off();
            $(this).on("click", function(e){
                remove_select();
                show_all_rows();
                set_processes_colors();
                $datacardTable.handsontable('getInstance').updateSettings({               
                    columns: colIDforDefault,
                    rowHeaders: rowHeaders,
                    colHeaders: bins
                });
                change_on_events_proc(true);
                change_on_events_bin(true);
                change_on_events_nuisances(true);
            });
        }
        else{
            $(this).off();
            var selectedProc = [];
            var binstoShow = [];
            binstoShow = getColDataIndicesOf(signals[procIndex], signals);
            for (var i = 0;i<binstoShow.length;i++){
                selectedProc.push(bins[binstoShow[i]["data"]]);
            }
            $(this).on("click", function(e){
                set_processes_colors(binstoShow);
                $datacardTable.handsontable('getInstance').updateSettings({              
                    columns: binstoShow,
                    colHeaders: selectedProc
                });
                change_on_events_proc(false);
                change_on_events_bin_proc(true, binstoShow, selectedProc);
                change_on_events_nuisances(false);
                hide_rows(binstoShow);
            });
        } 
        procIndex++;
    });
}

/*
   Function: change_on_events_nuisances
   
   Sorts datacard's handsontable by process.
      
   Parameters:

      firstClick - Boolean variable to identify if handsontable cell was clicked for the first time.
	  selected - jQuery element of selected nuisance.
*/
function change_on_events_nuisances(firstClick, selected){
    var rowIndex = 0;
    $("th").each(function() {
        if($(this).html().indexOf('<div class="relative">') < 0
        && $(this).html() != "process"){
            if(!firstClick || typeof(selected)!=='undefined'&& selected.children(0).html()==$(this).children(0).html()){
                $(this).off();
                $(this).on("click", function(e){
                    $(this).find('button:first').popover('destroy');
                    $(this).removeClass("selected");
                    show_all_rows();
                    set_processes_colors();
                    $datacardTable.handsontable('getInstance').updateSettings({               
                        columns: colIDforDefault,
                        rowHeaders: rowHeaders,
                        colHeaders: bins
                    });
                    change_on_events_nuisances(true);
                    change_on_events_bin(true);
                    change_on_events_proc(true);
                });
            }
            else{
                $(this).off();
                var cols = colIDtoShow[rowIndex];
                var colHead = colHeaderstoShow[rowIndex];
                $(this).on("click", function(e){
                    remove_select();
                    $(this).addClass("selected");
                    show_all_rows();
                    set_processes_colors(cols);
                    $datacardTable.handsontable('getInstance').updateSettings({              
                        columns: cols,
                        colHeaders: colHead
                    });
                    change_on_events_nuisances(true, $(this));
                    change_on_events_bin(false);
                    change_on_events_proc(false);
                });
            }
            add_menu_events($(this).find('button:first')); 
            rowIndex++;
        }
    });
}

/*
   Function: remove_select
   Remove selected style from datacard's handsontable when using search.
*/
function remove_select(){
    $("th").each(function() {
        if($(this).html().indexOf('<div class="relative">') < 0
        && $(this).html() != "process"){
            $(this).removeClass("selected");
        }   
    });
    destroy_popovers();
}

/*
   Function: destroy_popovers
   Remove popover element from datacard.
*/
function destroy_popovers(){
    var $popover = $("button.fa-bars").popover({
            selector: '[data-original-title=]'
    });
    $popover.popover('hide');
}

/*
   Function: add_menu_events
   
   Adds menu event on button click.
      
   Parameters:

      button - Selected jQuery button object.
*/
function add_menu_events(button){
    button.off();
    button.on("click", function(e){
        button.popover("destroy");
        var nuis = button.parent().children(0).first().html();
        var rowIndex = 0;
        $("th").each(function() {
            if($(this).html().indexOf('<div class="relative">') < 0
            && $(this).html() != "process"){
                if($(this).children(0).first().html() == nuis){
                    build_menu(button, rowIndex);
                    return false;
                }
                rowIndex++;
            }
        });
        e.stopPropagation();
    });
}

/*
   Function: count_header_width
   
   Searches for the longest jQuery th element and sets datacard's handsontable nuisance width.
      
   Parameters:

      th - Selected jQuery th object.
*/
function count_header_width(th){
    var size = th.children().first().textWidth();
    if (rowHeadersMaxWidth < size)
        rowHeadersMaxWidth = size;
}

/*
   Function: hide_rows
   
   Hides rows of nuisances without values.
      
   Parameters:

      cols - Object of datacard's handsontable columns shown.
*/
function hide_rows(cols){
    var handsonRows = $("table.htCore tbody tr");
    var toHide = [];
    for (var nu = 0;nu<nuisances.length;nu++){
        var count = 0;
        for (var i = 0; i<cols.length; i++){      
            if (nuisances[nu][cols[i]["data"]+2] == "-")
                count++;
        }
        if(count==cols.length)
            toHide.push(nu);
    }
    for (var i = 1; i < handsonRows.length; i++) {
        if (toHide.indexOf(i-1)>=0) {
            handsonRows[i].style.display = 'none';
        }
    }
}

/*
   Function: show_all_rows
   Shows all nuisances rows in datacard's handsontable.
*/
function show_all_rows(){
    $("table.htCore tbody tr").removeAttr("style");
}

/*
   Function: add_search_to_table
   
   Adds searching to datacard's handsontable.
      
   Parameters:

      tableData - Object of datacard's handsontable cell values.
*/
function add_search_to_table(tableData){
    $('#srch-term').off();
    $('#srch-term').on('keyup', function (event) {
        var value = ('' + this.value).toLowerCase();
        var td;

        if (value) {
            for (var row = 0; row < tableData.length; row++)
                for (var col = 0; col < tableData[row].length; col++) {
                    td = $datacardTable.handsontable('getCell', row, col);
                    if (('' + tableData[row][col]).toLowerCase().indexOf(value) > -1) 
                        td.className = 'pass'
                    else if (td != null)
                        td.className = '';
                }
        }else
            for (var row = 0; row < tableData.length; row++)
                for (var col = 0; col < tableData[row].length; col++) {
                    td = $datacardTable.handsontable('getCell', row, col);
                    if (td != null)
                        td.className = '';
                }
    });
}

/*
   Function: build_menu
   
   Creates popover menu from a nuisance button with sorting button choices.
      
   Parameters:

      button - Selected jQuery button object.
	  rowIndex - Selected nuisance index.
*/
function build_menu(button, rowIndex){
    button.off();
    button.on("click", function(e){
        e.stopPropagation();
    });
    button.popover({
        title: '<div class="text-center">Sort data</div>',
        content: function(){
            var menuIndex = 0;
            var html = "";
            html += '<div class="btn-group-vertical">';
            $.each(['By nuisance', 'By value', 'By process', 'By bin'/*, 'Others'*/], function (i, type) {
                html +=  '<button class="btn btn-default btn-xs" onclick="menu_sort('+menuIndex+','+rowIndex+');">'+type+'</button>';
                menuIndex++;
            });
            html += "</ul>";
            return html;
        },
        container: "body",
        placement: "top",
        html: true
    });
    button.popover('toggle');
}

/*
   Function: menu_sort
   
   On menu click opens Bootstrap Dialog with selected nuisance sorting.
      
   Parameters:

      menuIndex - Sorting menu id to show different table.
	  rowIndex - Selected nuisance index.
*/
function menu_sort(menuIndex, rowIndex){
    switch(menuIndex)
    {
    case 0:
        var cols = colIDtoShow[rowIndex];
        var content = colValuestoShow[rowIndex];
        var html = "";
        html += '<table class="table table-bordered text-left sortable">';
        html += "<thead><tr>";
        html += '<th data-defaultsort="asc" data-sortcolumn="0" data-sortkey="0-1">value</th>';
        html += '<th data-sortcolumn="1" data-sortkey="1-1">bin</th>';
        html += '<th data-sortcolumn="2" data-sortkey="2-1">process</th>';
        html += "</tr></thead><tbody>";
        for (var j=0;j<content.length;j++ ){
            html += "<tr>";
            if (content[j]+"".indexOf("<button") < 0)
                html += '<td data-value ="'+content[j]+'">'+content[j]+"</td>";
            else{
                html += '<td data-value ="'+j+'">'+content[j]+"</td>";
            }
            html += '<td data-value ="'+bins[cols[j].data]+'">'+bins[cols[j].data]+"</td>";
            html += '<td style="background-color:'+colors[cols[j].data]+'" data-value ="'+processes[cols[j].data]+'">'+processes[cols[j].data]+"</td>";
            html += "</tr>";
        }
        html += "</tbody></table>";

        BootstrapDialog.show({
            title: 'Datacard table sort by '+nuisances[rowIndex][0]+": "+nuisances[rowIndex][1]+' nuisance',
            message: function(){
                var $dialogTable = $(html);
                $dialogTable.find(":button").each(function(index, button) {
                    add_cell_dialog_event($(button));
                });
                return $dialogTable;
            }
        });
        break;
    case 1:
        var cols = colIDtoShow[rowIndex];
        var content = colValuestoShow[rowIndex];
        var othersValues = [];
        var othersIDs = [];
        for (var j = 0;j<content.length;j++){
            if(othersValues.indexOf(content[j]) == -1)
                othersValues.push(content[j]);
        }
        for (var j = 0;j<othersValues.length;j++){
            othersIDs.push(getColIndicesOf(othersValues[j], content));
        }
        var html = "";
        html += '<table class="table table-bordered text-left sortable">';
        html += "<thead><tr>";
        html += '<th data-defaultsort="asc" data-sortcolumn="0" data-sortkey="0-1">value</th>';
        html += '<th data-sortcolumn="1" data-sortkey="1-1">others (bin/process)</th>';
        html += "</tr></thead><tbody>";
        for (var j = 0;j<othersValues.length;j++){
            html += "<tr>";
            if (othersValues[j]+"".indexOf("<button") < 0)
                html += '<td data-value ="'+othersValues[j]+'">'+othersValues[j]+"</td>";
            else{
                html += '<td data-value ="'+j+'">'+othersValues[j]+"</td>";
            }
            var temp = "";
            for(var i = 0; i<othersIDs[j].length; i++){
                temp += bins[cols[othersIDs[j][i]].data]+'/' + processes[cols[othersIDs[j][i]].data]+'; ';
            }
            html += '<td data-value ="'+ temp +'">'+ temp+"</td>";
            html += "</tr>";
        }
        html += "</tbody></table>";

        BootstrapDialog.show({
            title: 'Datacard table sort by '+nuisances[rowIndex][0]+": "+nuisances[rowIndex][1]+' nuisance,<br>grouped by value',
            message: function(){
                var $dialogTable = $(html);
                $dialogTable.find(":button").each(function(index, button) {
                    add_cell_dialog_event($(button));
                });
                return $dialogTable;
            }
        });
        break;
    case 2:
        var cols = colIDtoShow[rowIndex];
        var content = colValuestoShow[rowIndex];
        var othersValues = [];
        var othersIDs = [];
        for (var j = 0;j<cols.length;j++){
            if(othersValues.indexOf(processes[cols[j].data]) == -1)
                othersValues.push(processes[cols[j].data]);
        }
        var proclist = [];
        for (var j = 0;j<cols.length;j++){
            proclist.push(processes[cols[j].data]);
        }
        for (var j = 0;j<othersValues.length;j++){
            othersIDs.push(getColIndicesOf(othersValues[j], proclist));
        }
        var html = "";
        html += '<table class="table table-bordered text-left sortable">';
        html += "<thead><tr>";
        html += '<th data-defaultsort="asc" data-sortcolumn="0" data-sortkey="0-1">process</th>';
        html += '<th data-sortcolumn="1" data-sortkey="1-1">others (bin/value)</th>';
        html += "</tr></thead><tbody>";
        for (var j = 0;j<othersValues.length;j++){
            html += "<tr>";
            html += '<td style="background-color:'+colors[processes.indexOf(othersValues[j])]+'" data-value ="'+othersValues[j]+'">'+othersValues[j]+"</td>";
            var temp = "";
            for(var i = 0; i<othersIDs[j].length; i++){
                temp += bins[cols[othersIDs[j][i]].data]+'/' + content[othersIDs[j][i]]+'; ';
            }
            html += '<td data-value ="'+ temp +'">'+ temp+"</td>";
            html += "</tr>";
        }
        html += "</tbody></table>";

        BootstrapDialog.show({
            title: 'Datacard table sort by '+nuisances[rowIndex][0]+": "+nuisances[rowIndex][1]+' nuisance,<br>grouped by process',
            message: function(){
                var $dialogTable = $(html);
                $dialogTable.find(":button").each(function(index, button) {
                    add_cell_dialog_event($(button));
                });
                return $dialogTable;
            }
        });
        break;
    case 3:
        var cols = colIDtoShow[rowIndex];
        var content = colValuestoShow[rowIndex];
        var othersValues = [];
        var othersIDs = [];
        for (var j = 0;j<cols.length;j++){
            if(othersValues.indexOf(bins[cols[j].data]) == -1)
                othersValues.push(bins[cols[j].data]);
        }
        var binlist = [];
        for (var j = 0;j<cols.length;j++){
            binlist.push(bins[cols[j].data]);
        }
        for (var j = 0;j<othersValues.length;j++){
            othersIDs.push(getColIndicesOf(othersValues[j], binlist));
        }
        var html = "";
        html += '<table class="table table-bordered text-left sortable">';
        html += "<thead><tr>";
        html += '<th data-defaultsort="asc" data-sortcolumn="0" data-sortkey="0-1">bin</th>';
        html += '<th data-sortcolumn="1" data-sortkey="1-1">others (process/value)</th>';
        html += "</tr></thead><tbody>";
        for (var j = 0;j<othersValues.length;j++){
            html += "<tr>";
            html += '<td data-value ="'+othersValues[j]+'">'+othersValues[j]+"</td>";
            var temp = "";
            for(var i = 0; i<othersIDs[j].length; i++){
                temp += processes[cols[othersIDs[j][i]].data]+'/' + content[othersIDs[j][i]]+'; ';
            }
            html += '<td data-value ="'+ temp +'">'+ temp+"</td>";
            html += "</tr>";
        }
        html += "</tbody></table>";

        BootstrapDialog.show({
            title: 'Datacard table sort by '+nuisances[rowIndex][0]+": "+nuisances[rowIndex][1]+' nuisance,<br>grouped by bin',
            message: function(){
                var $dialogTable = $(html);
                $dialogTable.find(":button").each(function(index, button) {
                    add_cell_dialog_event($(button));
                });
                return $dialogTable;
            }
        });
        break;
    case 4:
        //TODO OTHERS
        break;
    default:
        //code to be executed if n is different from case 1 and 2
    }
}