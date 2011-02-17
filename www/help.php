<?php
if (true) {
    ini_set('display_errors', 'On');
    error_reporting(E_ALL);
}
$relDir = ".";
include_once("$relDir/Menu.php");
include_once("$relDir/Page.php");

$menu = new Menu();
$page = new Page("LSST Pipetest", "LSST Pipe Test Help", $menu);


# Now write the main help content page.
include_once("Html.php");

#######################################
# Steps to add new data as <ol>
$ulNewData = new OrderedList();
$ulNewData->addItem("Add a fooTestDataNNN directory to a directory visible in your TESTBED_PATH.");
$ulNewData->addItem("Put your data: bias/flat/raw/registry into the new direcortory.");
$ulNewData->addItem("Get checksums for your data directory. (see bin/writeTestDataManifest.py).");
$ulNewData->addItem("Add a class 'FooTestData' which inherits from TestData.");
# class you include and error message to explain origin of data

$page->appendContent("<h2>How To Add New Data</h2>\n");
$page->appendContent($ulNewData->write());


#######################################
# Steps to add new test as <ul>
$ulNewTest = new OrderedList();
$ulNewTest->addItem("Add any new data as described above, but try to use existing data if possible.");
$ulNewTest->addItem("Follow tests/psfPhotometry.py as an example.");

$page->appendContent("<h2>How To Add a New Test</h2>\n");
$page->appendContent($ulNewTest->write());

########################################
# write the page
echo $page;
