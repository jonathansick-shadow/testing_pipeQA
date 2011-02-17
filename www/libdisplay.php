<?php

include_once("Html.php");
include_once("libdb.php");


######################################
# true/false color
# color text green if true, else red 
######################################
function tfColor($string, $tf) {
    if ( $tf ) {
        $colorout = "<font color=\"#00aa00\">$string</font>";
    } else {
        $colorout = "<font color=\"#880000\">$string</font>";
    }
    return $colorout;
}


#######################################
#
#
#######################################
function getCurrentUriDir() {
    $uri = $_SERVER['REQUEST_URI'];
    $dirs = preg_split("/\//", $uri);
    $dir = $dirs[count($dirs)-2];
    return $dir;
}


function verifyTest($value, $lo, $hi) {
    $pass = true;  # default true (ie. no limits were set)
    if ($lo and $hi) {
	$pass = ($value >= $lo and $value <=$hi);
    }
    if ($lo and !$hi) {
	$pass = ($value >= $lo);
    }
    if (!$lo and $hi) {
	$pass = ($value <= $hi);
    }
    return $pass;
}



function writeTable_ListOfTestResults($testDir) {

    $table = new Table("width=\"90%\"");
    
    $table->addHeader(array("Label", "Timestamp", "LowerLimit", "Value", "UpperLimit", "Comment"));

    global $dbFile;
    $db = connect();
    $cmd = "select * from summary";
    $result = $db->query($cmd);
    
    foreach ($result as $r) {
	list($test, $lo, $value, $hi, $comment) =
	    array($r['label'], $r['lowerlimit'], $r['value'], $r['upperlimit'], $r['comment']);
	
	$pass = verifyTest($value, $lo, $hi);
	if (!$lo) { $lo = "None"; }
	if (!$hi) { $hi = "None"; }

	if (!$pass) {
	    $test .= " <a href=\"backtrace.php?label=$test\">Backtrace</a>";
	}
	$mtime = date("Y-m_d H:i:s", $r['entrytime']);
	
	$table->addRow(array($test, $mtime, $lo, tfColor($value, $pass), $hi, $comment));
    }
    $db = NULL;
    return $table->write();
    
}
function displayTable_ListOfTestResults($testDir) {
    echo writeTable_ListOfTestResults($testDir);
}



function writeTable_OneTestResult($testDir, $label) {

    if (empty($label)) {
	return "<h2>No test label specified. Cannot display test result.</h2><br/>\n";
    }
    
    $table = new Table("width=\"90%\"");
    
    $table->addHeader(array("Label", "Timestamp", "LowerLimit", "Value", "UpperLimit", "Comment"));

    global $dbFile;
    #$mtime = date("Y-m_d H:i:s", filemtime("$testDir/$dbFile"));
    $db = connect();
    $cmd = "select * from summary where label='$label'";
    $result = $db->query($cmd);
    
    foreach ($result as $r) {
	list($test, $timestamp, $lo, $value, $hi, $comment, $backtrace) =
	    array($r['label'], $r['entrytime'], $r['lowerlimit'], $r['value'], $r['upperlimit'],
		  $r['comment'], $r['backtrace']);
	
	$pass = verifyTest($value, $lo, $hi);
	if (!$lo) { $lo = "None"; }
	if (!$hi) { $hi = "None"; }

	$table->addRow(array($test, date("Y-m-d H:i:s", $timestamp),
			     $lo, tfColor($value, $pass), $hi, $comment));
    }
    $db = NULL;

    return $table->write();
}
function write_OneBacktrace($testDir, $label) {

    $out = "<h2>Backtrace</h2><br/>\n";
    if (empty($label)) {
	return "<b>No test label specified. Cannot display backtrace.</b><br/>\n";
    }
    
    global $dbFile;
    $db = connect();
    $cmd = "select * from summary where label='$label'";
    $result = $db->query($cmd);

    $backtrace = "";
    foreach ($result as $r) {
	$backtrace .= $r['backtrace'];
    }
    $db = NULL;

    $out .= preg_replace("/\n/", "<br/>\n", $backtrace);
    $out = preg_replace("/(\t|\s{4})/", "&nbsp;&nbsp;", $out);
    
    return $out;
    
}
function displayTable_OneTestResult($testDir, $label) {
    echo writeTable_OneTestResult($testDir);
}






function writeFigures($testDir) {
    
    $out = "";
    $d = @dir("$testDir");
    while(false !== ($f = $d->read())) {
    	if (! ereg(".(png|PNG|jpg|JPG)", $f)) { continue; }

	# get the image path
    	$path = "$testDir/$f";
    	$mtime = date("Y-m_d H:i:s", filemtime($path));

	# get the caption
	$db = connect();
	$cmd = "select caption from figure where filename = '$f'";
	$result = $db->query($cmd)->fetchColumn();

	$img = new Table();
	$img->addRow(array("<img src=\"$path\">"));
	$img->addRow(array($result));
	$img->addRow(array("Timestamp: $mtime"));
	$out .= $img->write();
    }
    return $out;
}
function displayFigures($testDir) {
    echo writeFigures($testDir);
}




function summarizeTest($testDir) {
    $summary = array();

    global $dbFile;
    #$mtime = date("Y-m_d H:i:s", filemtime("$testDir/$dbFile"));

    if (!file_exists($dbFile)){
	return "";
    }
    $db = connect($testDir);
    $testCmd = "select count(*) from summary";
    $nTest = $db->query($testCmd)->fetchColumn();
    
    $passCmd = "select * from summary";
    $results = $db->query($passCmd);
    $nPass = 0;
    $timestamp = 0;
    foreach($results as $result) {
	if (verifyTest($result['value'], $result['lowerlimit'], $result['upperlimit'])){
	    $nPass += 1;
	}
	$timestamp = $result['entrytime'];
    }

    $ret = array();
    $ret['name'] = $testDir;
    $ret['entrytime'] = $timestamp;
    $ret['ntest'] = $nTest;
    $ret['npass'] = $nPass;
    return $ret;
}


function writeTable_SummarizeAllTests() {
    $dir = "./";
    
    ## go through all directories and look for .summary files
    $d = @dir($dir) or dir("");

    $table = new Table("width=\"90%\"");
    $table->addHeader(array("Test", "mtime", "No. Tests", "No. Passed"));
    while(false !== ($testDir = $d->read())) { 
	if ( ereg("^\.", $testDir) or ! is_dir("$testDir")) {
	    continue;
	}
	$summ = summarizeTest($testDir);
	$testLink = "<a href=\"$testDir\">$testDir</a>";
	$passLink = tfColor($summ['npass'], ($summ['npass']==$summ['ntest']));
	$table->addRow(array($testLink, date("Y-m-d H:i:s", $summ['entrytime']), $summ['ntest'], $passLink));
    }
    return $table->write();
    
}

function displayTable_SummarizeAllTests() {
    echo writeTable_SummarizeAllTests();
}



function writeTable_Logs() {
    
    $db = connect();

    # first get the tables ... one for each ccd run
    $cmd = "select name from sqlite_sequence where name like 'log%'";
    $dbtables = $db->query($cmd);

    # make links at the top of the page
    $tables = "";
    $ul = new UnorderedList();
    foreach ($dbtables as $dbtable) {

	$name = $dbtable['name'];
	$ul->addItem("<a href=\"#$name\">$name</a>");
	
	$cmd = "select * from $name";
	$logs = $db->query($cmd);

	$tables .= "<h2 id=\"$name\">$name</h2><br/>";
	
	$table = new Table("width=\"80%\"");
	$table->addHeader(array("Module", "Message", "Date", "Level"));
	foreach ($logs as $log) {

	    # check for tracebacks from TestData
	    $module = $log['module'];
	    $msg = $log['message'];
	    if (ereg("testQA.TestData$", $module)) {
		# get the idString from the message
		$idString = preg_replace("/:.*/", "", $msg);
		$module .= " <a href=\"backtrace.php?label=$idString\">Backtrace</a>";
	    }
	    $table->addRow(array($module, $msg, $log['date'], $log['level']));
	}
	$tables .= $table->write();
    }
    $contents = "<h2>Data Used in This Test</h2><br/>" . $ul->write() . "<br/><br/>";
    return $contents . $tables;
}

function displayTable_Logs() {
    echo writeTable_Logs();
}


function writeTable_EupsSetups() {

    $db = connect();

    # first get the tables ... one for each ccd run
    $cmd = "select name from sqlite_sequence where name like 'eups%'";
    $dbtables = $db->query($cmd);

    # make links at the top of the page
    $tables = "";
    $ul = new UnorderedList();
    foreach ($dbtables as $dbtable) {

	$name = $dbtable['name'];
	$ul->addItem("<a href=\"#$name\">$name</a>");
	
	$cmd = "select * from $name";
	$logs = $db->query($cmd);

	$tables .= "<h2 id=\"$name\">$name</h2><br/>";
	
	$table = new Table("width=\"80%\"");
	$table->addHeader(array("Product", "Version", "Timestamp"));
	foreach ($logs as $log) {
	    $table->addRow(array($log['product'],$log['version'],date("Y-m-d H:i:s", $log['entrytime'])));
	}
	$tables .= $table->write();
    }
    $contents = "<h2>Data Sets Used in This Test</h2><br/>" . $ul->write() . "<br/><br/>";
    return $contents . $tables;

}

function displayTable_EupsSetups() {
    echo writeTable_EupsSetups();
}



?>