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
$page = new Page("LSST Pipetest", "EUPS: setup products", $menu);

$page->appendContent(writeTable_EupsSetups());

echo $page;
