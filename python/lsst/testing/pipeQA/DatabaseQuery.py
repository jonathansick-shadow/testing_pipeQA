import MySQLdb

# Base class
class DatabaseInterface():
    def __init__(self):
        pass

# LSST specific interface
class LsstSimDbInterface(DatabaseInterface):
    # Mapping from filter names to database names
    filterMap = { "u" : 0, "g" : 1, "r" : 2, "i" : 3, "z" : 4 }
    mySqlHost   = 'lsst10.ncsa.uiuc.edu'
    mySqlUser   = 'abecker'
    mySqlPasswd = 'Timp@ne80x'

    def __init__(self, mySqlDb):
        DatabaseInterface.__init__(self)
        
        self.db     = MySQLdb.connect(host=LsstSimDbInterface.mySqlHost,
                                      user=LsstSimDbInterface.mySqlUser,
                                      passwd=LsstSimDbInterface.mySqlPasswd,
                                      db=mySqlDb)
        self.cursor = self.db.cursor()

    def execute(self, sql, verbose = False):
        if verbose:
            print sql
        self.cursor.execute(sql)
        return self.cursor.fetchall()

    
