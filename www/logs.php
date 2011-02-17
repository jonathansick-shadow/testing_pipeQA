<?php
if (true) {
    ini_set('display_errors', 'On');
    error_reporting(E_ALL);
}
$relDir = "..";
include_once("$relDir/Menu.php");
include_once("$relDir/Page.php");
include_once("$relDir/libdisplay.php");

$menu = new TestMenu();
$page = new Page("LSST Pipetest", "Test Logs", $menu);
$page->appendContent(writeTable_Logs());

echo $page;
