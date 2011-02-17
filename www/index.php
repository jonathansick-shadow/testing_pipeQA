<?php
if (true) {
    ini_set('display_errors', 'On');
    error_reporting(E_ALL);
}
$relDir = ".";
include("$relDir/Menu.php");
include("$relDir/Page.php");
include("$relDir/libdisplay.php");

$menu = new Menu();
$page = new Page("LSST Pipetest", "LSST Pipe Test Summary", $menu);
$page->appendContent(writeTable_SummarizeAllTests());

echo $page;
