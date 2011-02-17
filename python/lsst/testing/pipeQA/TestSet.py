import sys
import traceback
import os
import re
import inspect
import sqlite3
import stat
import eups

from lsst.testing.pipeQA.LogConverter  import LogFileConverter
from lsst.testing.pipeQA.TestFailError import TestFailError

class TestSet(object):
    
    def __init__(self, mainDisplayFull=False):
        """Constructor to create a TestSet object for a new suite of tests."""
        
        wwwBase = "www"
        testfileName = inspect.stack()[-1][1]
        self.testfileBase = re.sub(".py", "", os.path.split(testfileName)[1])
        self.wwwDir = os.path.join(wwwBase, self.testfileBase)

        if not os.path.exists(self.wwwDir):
            os.mkdir(self.wwwDir)


        # create symlinks to the test-specific pages
        toLink = [
            ["summary.php",   "summary.php"],
            ["logs.php",      "logs.php"],
            ["sdqa.php",      "sdqa.php"],
            ["eups.php",      "eups.php"],
            ["redirect.php",  "index.php"],
            ["backtrace.php", "backtrace.php"],
            ]
        for pair in toLink:
            srcFile = os.path.join("..", pair[0])
            symlink = os.path.join(self.wwwDir, pair[1])
            if not os.path.exists(symlink):
                os.symlink(srcFile, symlink)

                
        # connect to the db and create the tables
        self.dbFile = os.path.join(self.wwwDir, "db.sqlite3")
        self.conn = sqlite3.connect(self.dbFile)
        self.curs = self.conn.cursor()
        self.summTable, self.figTable, self.eupsTable = "summary", "figure", "eups"
        self.tables = {
            self.summTable : ["label text unique", "value double",
                              "lowerlimit double", "upperlimit double", "comment text",
                              "backtrace text"],
            self.figTable  : ["filename text", "caption text"],
            }

        self.stdKeys = ["id integer primary key autoincrement", "entrytime timestamp"]
        for k, v in self.tables.items():
            keys = self.stdKeys + v
            cmd = "create table if not exists " + k + " ("+",".join(keys)+")"
            self.conn.execute(cmd)

        
        self.conn.commit()

        self.tests = []
        
    def __del__(self):
        self.conn.close()


        
    def _insertOrUpdate(self, table, replacements, selectKeys):
        """Insert entries into a database table, overwrite if they already exist."""
        
        # there must be a better sql way to do this ... but my sql-foo is weak
        # we want to overwrite entries if they exist, or insert them if they don't
        
        # delete the rows which match the selectKeys
        where = []
        for key in selectKeys:
            if isinstance(replacements[key], str):
                where.append(key + "='" + replacements[key] + "'")
            else:
                where.append(key + "=" + str(replacements[key]))
        where = " where " + " and ".join(where)
        
        cmd = "delete from " + table + " " + where
        self.curs.execute(cmd)

        
        # insert the new data
        keys = []
        values = []
        for k,v in replacements.items():
            keys.append(k)
            values.append(v)
        values = tuple(values)
        inlist = " (id, entrytime,"+ ",".join(keys) + ") "
        qmark = " (NULL, strftime('%s', 'now')," + ",".join("?"*len(values)) + ")"
        cmd = "insert into "+table+inlist + " values " + qmark

        self.curs.execute(cmd, values)
        self.conn.commit()

    
    def addTest(self, label, value, limits, comment):
        """Add a test to this testing suite."""
        
        self.tests.append(label)

        # grab a traceback for failed tests
        backtrace = ""
        try:
            if (value < limits[0] or value > limits[1]):
                raise TestFailError("Failed test '"+label+"': " +
                                    "value '" + str(value) + "' not in range '" + str(limits)+"'.")
        except TestFailError, e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            backtrace = "".join(traceback.format_stack()[:-1]) + "\n" + str(e)
            
        # enter the test in the db
        keys = [x.split()[0] for x in self.tables[self.summTable]]
        replacements = dict( zip(keys, [label, value, limits[0], limits[1], comment, backtrace]) )
        self._insertOrUpdate(self.summTable, replacements, ['label'])

        
    def importExceptionDict(self, exceptDict):
        """Given a dictionary of exceptions from TestData object, add the entries to the db."""
        keys = sorted(exceptDict.keys())
        for key in keys:
            tablekeys = [x.split()[0] for x in self.tables[self.summTable]]
            replacements = dict( zip(tablekeys, [key, 0, 1, 1, "Uncaught exception", exceptDict[key]]) )
            self._insertOrUpdate(self.summTable, replacements, ['label'])

        
    def addFigure(self, fig, filename, caption):
        """Add a figure to this test suite."""
        path = os.path.join(self.wwwDir, filename)
        fig.savefig(path)
        
        keys = [x.split()[0] for x in self.tables[self.figTable]]
        replacements = dict( zip(keys, [filename, caption]))
        self._insertOrUpdate(self.figTable, replacements, ['filename'])
        
        
    def importLogs(self, logFiles):
        """Import logs from logFiles output by pipette."""
        
        # convert our ascii logfile to a sqlite3 database    
        def importLog(logFile):
            base = os.path.basename(logFile)
            table = "log_" + re.sub(".log", "", base)
            converter = LogFileConverter(logFile)
            converter.writeSqlite3Table(self.dbFile, table)

        # allow a list of filenames to be provided, or just a single filename
        if isinstance(logFiles, list):
            for logFile in logFiles:
                importLog(logFile)
        else:
            importLog(logFile)

            
    def importEupsSetups(self, eupsSetupFiles):
        """Import the EUPS setup packages from files written by TestData object during pipette run."""

        # note that this only works if we ran the test ourselves.

        def importEups(eupsFile):
            base = os.path.basename(eupsFile)
            table = "eups_" + re.sub(".eups", "", base)
            setups = []
            fp = open(eupsFile, 'r')
            for line in fp.readlines():
                setups.append(line.split())
            fp.close()

            mykeys = ["product text", "version text"]
            keys = self.stdKeys + mykeys

            cmd = "create table if not exists " + table + " ("+",".join(keys)+")"
            self.curs.execute(cmd)
            self.conn.commit()

            for setup in setups:
                product, version = setup
                mykeys = [x.split()[0] for x in mykeys]
                replacements = dict( zip(mykeys, [product, version]))
                self._insertOrUpdate(table, replacements, ['product'])
            
        # allow a list of files or just one
        if isinstance(eupsSetupFiles, list):
            for eupsSetupFile in eupsSetupFiles:
                importEups(eupsSetupFile)
        else:
            importEups(eupsSetupFile)
                
