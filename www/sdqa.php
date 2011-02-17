<?php
if (true) {
    ini_set('display_errors', 'On');
    error_reporting(E_ALL);
}
$relDir = "..";
include("$relDir/Menu.php");
include("$relDir/Page.php");

$menu = new TestMenu();
$page = new Page("LSST Pipetest", "SDQA Information", $menu);
$page->appendContent("Not Yet Implemented.\n");

echo $page;
