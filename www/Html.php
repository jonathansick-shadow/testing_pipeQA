<?php

class Div {
    private $_tagAttribs;
    private $_content = array();
    
    public function __construct($tagAttribs) {
	$this->_tagAttribs = $tagAttribs;
    }
    public function append($content) {
	$this->_content[] = $content;
    }
    public function __toString() {
	$div = "<div $this->_tagAttribs>\n";
	foreach ($this->_content as $content) {
	    $div .= $content;
	}
	$div .= "</div>\n";
	return $div;
    }
    public function write() { return $this->__toString(); }
}

class HtmlList {

    private $_tagAttribs;
    private $_type;
    private $_items = array();
    
    public function __construct($type = "ul", $tagAttribs = "") {
	$this->_tagAttribs = $tagAttribs;
	$this->_type = $type;
    }
    
    public function addItem($item, $tagAttribs = "") {
	$this->_items[] = array($item, $tagAttribs);
    }
    
    public function __toString() {
	$s = "<".$this->_type." ".$this->_tagAttribs.">\n";
	foreach ($this->_items as $item) {
	    $it = $item[0];
	    $tagAttrib = $item[1];
	    $s .= "<li ".$tagAttrib .">$it</li>\n";
	}
	$s .= "</".$this->_type.">\n";
	return $s;
    }
    public function write() { return $this->__toString(); }
}

class OrderedList extends HtmlList {
    public function __construct($tagAttribs="") {
	parent::__construct("ol", $tagAttribs);
    }
}
class UnorderedList extends HtmlList {
    public function __construct($tagAttribs="") {
	parent::__construct("ul", $tagAttribs);
    }
}


class Table {
    private $_header = array();
    private $_rows = array();
    private $_tagAttribs;
    
    public function __construct($tagAttribs = "") {
	$this->_tagAttribs = $tagAttribs;
    }
    public function addHeader($array) {
	$this->_header = $array;
    }
    public function addRow($array) {
	$this->_rows[] = $array;
    }
    public function __toString() {
	$s = "<table ".$this->_tagAttribs.">\n";
	$s .= $this->_writeRow($this->_header);
	foreach ($this->_rows as $row) {
	    $s .= $this->_writeRow($row);
	}
	$s .= "</table>\n";
	return $s;
    }
    public function write() { return $this->__toString(); }
    private function _writeRow($array, $td="td") {
	$row = "<tr>";
	foreach ($array as $entry) {
	    $row .= "<$td>$entry</$td>";
	}
	$row .= "</tr>\n";
	return $row;
    }
}
