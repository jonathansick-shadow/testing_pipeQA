import os
import MySQLdb
import lsst.pex.policy as pexPolicy
import time
from lsst.pex.logging import Trace

class DatabaseIdentity:
    """
    Requires file that looks like:
    
    cat ~/.lsst/db-auth.paf
    database: {
        authInfo: {
            host: lsst10.ncsa.uiuc.edu
            user: "XXXX"
            password: "YYYY"
        }
    }
    """
    def __init__(self, mySqlDb):
        self.mySqlDb   = mySqlDb
        self.loadId()

    def loadId(self):
        dbAuth = os.path.join(os.environ["HOME"], ".lsst", "db-auth.paf")
        policy = pexPolicy.Policy(dbAuth)
        authPolicy = policy.get("database").get("authInfo")
        self.mySqlUser = authPolicy.get("user")
        self.mySqlHost = authPolicy.get("host")
        self.mySqlPasswd = authPolicy.get("password")
        

# Base class
class DatabaseInterface():
    def __init__(self):
        pass


# LSST specific interface
class LsstSimDbInterface(DatabaseInterface):
    # Mapping from filter names to database names
    filterMap = { "u" : 0, "g" : 1, "r" : 2, "i" : 3, "z" : 4 }

    def __init__(self, dbId):
        """
        @param dbId  A databaseIdentity object contain connection information
        """
        DatabaseInterface.__init__(self)
        
        self.db     = MySQLdb.connect(host = dbId.mySqlHost,
                                      db = dbId.mySqlDb,
                                      user = dbId.mySqlUser,
                                      passwd = dbId.mySqlPasswd)
        self.cursor = self.db.cursor()

    def execute(self, sql):
        """Execute an sql command

        @param sql Command to be executed.
        """

        Trace("lsst.testing.pipeQA.LsstSimDbInterface", 3, "Executing: %s" % (sql))
        t0 = time.time()
        self.cursor.execute(sql)
        results = self.cursor.fetchall()
        t1 = time.time()
        Trace("lsst.testing.pipeQA.LsstSimDbInterface", 2, "Time for SQL query: %.1f s" % (t1-t0))
        return results

    
