<?php

$dbFile = "db.sqlite3";

function connect($dir=".") {
    global $dbFile;
    try {
	$db = new PDO("sqlite:$dir/$dbFile");
    } catch(PDOException $e) {
	print 'Exception : '.$e->getMessage();
    }
    return $db;
}

?>