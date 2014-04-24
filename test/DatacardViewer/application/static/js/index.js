$(function() {
    $("#tabs").tabs();
    $("#tabs li").click(function() {
	    $('li').removeClass('active');
	    $(this).addClass("active");
	});
});
$(document).on("click", function (e) {
    var $popover = $("button.fa-bars").popover({
        selector: '[data-original-title=]'
    });
    var $target = $(e.target),
        isPopover = $(e.target).is('[data-original-title=]'),
        inPopover = $(e.target).closest('.popover').length > 0
	
    //hide only if clicked on button or inside popover
    if (!isPopover && !inPopover)
    	$popover.popover('hide');
});