
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
var settings_mass = 120.7;
var histograms = {};
var histogramsWidth = {};

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
}

function get_colors(signals){
    var backgrounds = ["#DDA0DD", "#FCF8E3", "#DFF0D8", "#D9EDF7", "#F2DEDE", "#FFE4C4",
                        "#DEB887", "#BDB76B", "#E9967A", "#8FBC8F", "#CD5C5C", "F08080",
                        "#20B2AA", "#778899", "#191970", "#FF4500", "#D8BFD8", "#708090"];
    var colors = [];
    for(var i=0;i<signals.length;i++)
        colors.push(backgrounds[signals[i]]);
    return colors;
}

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


function getColDataIndicesOf(searchArr, arr) {
    var indices = [];
    for (var i = 0; i< arr.length; i++){
        if (arr[i]==searchArr)
            indices.push({"data":i});
    }
    return indices;
}

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

//from loadJSRootIO.js
function init_JSRootIO(){
    assertPrerequisitesAndRead();
    //TODO mass settings
    //set_settings();
}

function show_datacard(data){
    init_variables();
    parse_datacard(data);
    generate_datacard();
}

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

var processRenderer = function (instance, td, row, col, prop, value, cellProperties) {
    Handsontable.renderers.TextRenderer.apply(this, arguments);
    $(td).css({
        background: colorsToShow[col]
    });
};
var descriptionRenderer = function (instance, td, row, col, prop, value, cellProperties) {
    var escaped = Handsontable.helper.stringify(value);
    escaped = strip_tags(escaped, '<em><b><a>'); 
    td.innerHTML = escaped;
    return td;
};

var cellButtonRenderer = function (instance, td, row, col, prop, value, cellProperties) {
    if((value+"").indexOf("<button")>=0){
        td.innerHTML = Handsontable.helper.stringify((value+"").insert(value.indexOf('>'), ' style="width:'+$(td).width()+'px"'));
        add_cell_dialog_event($(td).find('button:first'));
    }else
        td.innerHTML = value;
};
//todo
function add_cell_dialog_event(button){
    button.off();
    button.on("click", function(e){
        build_cell_dialog(button);
        e.stopPropagation();
    });
}

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
                    var html = "<b>Nominal: </b><span style='color:red;'>red</span>, ";
                    html+="<b>Up</b>: <span style='color:blue;'>blue</span>, ";
                    html+="<b>Down</b>: <span style='color:green;'>green</span>.";
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
     
function change_headers_width(){
    $(".handsontable col.rowHeader ").css( "width", rowHeadersMaxWidth + 50);
    $datacardTable.handsontable('getInstance').render();
}

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

function remove_select(){
    $("th").each(function() {
        if($(this).html().indexOf('<div class="relative">') < 0
        && $(this).html() != "process"){
            $(this).removeClass("selected");
        }   
    });
    destroy_popovers();
}

function destroy_popovers(){
    var $popover = $("button.fa-bars").popover({
            selector: '[data-original-title=]'
    });
    $popover.popover('hide');
}

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

function count_header_width(th){
    var size = th.children().first().textWidth();
    if (rowHeadersMaxWidth < size)
        rowHeadersMaxWidth = size;
}

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

function show_all_rows(){
    $("table.htCore tbody tr").removeAttr("style");
}

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
//todo buttons build_cell_dialog
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
            html += '<td data-value ="'+content[j]+'">'+content[j]+"</td>";
            html += '<td data-value ="'+bins[cols[j].data]+'">'+bins[cols[j].data]+"</td>";
            html += '<td style="background-color:'+colors[cols[j].data]+'" data-value ="'+processes[cols[j].data]+'">'+processes[cols[j].data]+"</td>";
            html += "</tr>";
        }
        html += "</tbody></table>";

        BootstrapDialog.show({
            title: 'Datacard table sort by '+nuisances[rowIndex][0]+": "+nuisances[rowIndex][1]+' nuisance',
            message: html
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
            html += '<td data-value ="'+othersValues[j]+'">'+othersValues[j]+"</td>";
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
            message: html
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
            message: html
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
            message: html
        });
        break;
    case 4:
        //execute code block 1
        break;
    default:
        //code to be executed if n is different from case 1 and 2
    }
}