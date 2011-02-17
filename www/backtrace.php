<?php
if (true) {
    ini_set('display_errors', 'On');
    error_reporting(E_ALL);
}
$relDir = "..";
include("$relDir/Menu.php");
include("$relDir/Page.php");
include("$relDir/libdisplay.php");

$menu = new SpecificTestMenu();
$page = new Page("LSST Pipetest", "LSST Pipe Test Summary", $menu);

$page->appendContent("<h2>".getCurrentUriDir()."</h2><br/>\n");
$page->appendContent(writeTable_OneTestResult(".", $_GET['label']) . "<br/>\n");
$page->appendContent(write_OneBackTrace(".", $_GET['label']));

echo $page;
