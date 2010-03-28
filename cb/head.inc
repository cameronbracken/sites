document.write('<link href="stc.css" rel="stylesheet" type="text/css" />            ')
document.write('<link rel="shortcut icon" href="images/favicon.ico" type="image/ico"> ')

function lunchboxOpen(lunchID) {
	var state = document.getElementById('lunch_' + lunchID).style.display;
	
	if (state == "block") { 
		document.getElementById('lunch_' + lunchID).style.display = "none";
	}
	else
	{
		document.getElementById('lunch_' + lunchID).style.display = "block";
	}
}
