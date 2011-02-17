<?php
if (true) {
    ini_set('display_errors', 'On');
    error_reporting(E_ALL);
}
$relDir = "..";
include("$relDir/Menu.php");
include("$relDir/Page.php");
include("$relDir/libdisplay.php");

$menu = new TestMenu();
$page = new Page("LSST Pipetest", "LSST Pipe Test Summary", $menu);

$page->appendContent("<h2>".getCurrentUriDir()."</h2>\n");
$page->appendContent(writeTable_ListOfTestResults("."));
$page->appendContent(writeFigures("."));

echo $page;
