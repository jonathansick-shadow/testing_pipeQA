import os, re, sqlite3

class LogFileConverter (object):
    """Convert a pex_logging ascii log file to another format."""
    
    def __init__(self, logFile):

        self.logFile = logFile

        path, filename = os.path.split(logFile)
        idString = re.sub(".log", "", path)

        # \n delimits fields
        # \n\n delimits records (ie. a blank line)

        self.allLogs = []
        thisLog = []
        fp = open(self.logFile, 'r')
        for line in fp.readlines():

            # split on first ':' with a regex match
            m = re.search("^([^:]+):(.*)$", line)

            # if it matches, determine if it's -> module: message
            #                             or   -> LABEL: value
            if m:
                x = m.groups()
                if re.search("^\s*[A-Z]+$", x[0]):
                    thisLog.append(x[1].strip())
                else:
                    thisLog.append(x[0].strip())
                    thisLog.append(x[1].strip())
                
            # non-match means end of record
            else:
                self.allLogs.append(thisLog)
                thisLog = []


        
    def writeSqlite3Table(self, dbFile, table):
        """Write the log file to a sqlite3 database table"""

        table = re.sub("-", "_", table)
        
        conn = sqlite3.connect(dbFile)

        # add the table
        keys = ["id integer primary key autoincrement", "module text", "message text",
                "date text", "timestamp text", "level text"]
        cmd = "create table if not exists " + table + " ("+",".join(keys)+")"
        conn.execute(cmd)

        # insert the data
        for thisLog in self.allLogs:
            cmd = "insert into " + table + " values (NULL, ?, ?, ?, ?, ?)"
            conn.execute(cmd, tuple(thisLog))
        conn.commit()
        
        conn.close()
            

