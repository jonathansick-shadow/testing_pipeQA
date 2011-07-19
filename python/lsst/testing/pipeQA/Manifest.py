import os, re
import datetime
from lsst.testing.pipeQA.Checksum import Checksum

# could inherit from Dict perhaps
class ManifestHeader(object):
    """The header metadata for a manifest file. """

    def __init__(self, strings=[]):
        self.entries = {}

        # define required keys
        self.requiredKeys = ['hashtype', 'date']
        for k in self.requiredKeys:
            self.entries[k] = None

        # if we were given a header string, parse it
        if strings:
            for line in strings:
                m = re.search("#\s+(\w+)\s+(.*)$", line)
                if m:
                    k, v = m.groups()
                    self.entries[k] = v
        else:
            pass

    def set(self, key, value):
        """Set a header key,value pair."""
        self.entries[key] = value

    def get(self, key):
        """Get a header entry."""
        return self.entries[key]

    
    def verify(self):
        """Verify required keys are present in header."""
        
        # check for missing keys
        missingKeys = []
        for k in self.requiredKeys:
            if self.entries[k] is None:
                missingKeys.append(k)
        return missingKeys

    
    def write(self):
        """Write header to a string."""
        
        missingKeys = self.verify()
        if missingKeys:
            print "Warning: required manifest header keys missing:\n" + "\n".join(missingKeys)

        outstring = ""
        for k,v in self.entries.items():
            outstring += "# " + k + " " + v + "\n"
        return outstring


    
class Manifest(object):
    """A class to manage data file locations and checksums."""

    def __init__(self, testdataDir):
        """
        @param testdataDir   directory containing data.
        """
        self.testdataDir  = testdataDir
        self.checksums    = {}
        self.filepaths    = []
        self.header       = ManifestHeader([])

        self.manifest = os.path.join(self.testdataDir, "manifest")
        self.haveManifest = False
        if os.path.exists(self.manifest):
            self.haveManifest = True

        
    def read(self):
        """Load a manifest file. """
        
        missingInputs = []
        failedMd5s = []

        if not self.haveManifest:
            print "Unable to read manifest.  File "+self.manifest+" does not exist."
            return

        fp = open(self.manifest, 'r')
        lines = fp.readlines()
        self.header = ManifestHeader(lines)

        self.checksums = {}
        self.filepaths = []
        for line in lines:
            if re.search("^#", line):
                continue
            if not re.search("[^\s]+\s+[^\s]+", line):
                continue
            filename, checksum = line.split()
            filepath = os.path.join(self.testdataDir, filename)
            self.checksums[filepath] = checksum
            self.filepaths.append(filepath)

        fp.close()

        
    def getHeader(self, key):
        """Get a specified key from the manifest header."""
        
        if self.haveManifest:
            return self.header.get(key)
        else:
            return None

    def verifyExists(self):
        """Verify that files we're managing exist. """
        
        missingInputs = []
        for filepath in self.filepaths:
            if not os.path.exists(filepath):
                missingInputs.append(filepath)
        return missingInputs

    def verifyChecksum(self):
        """Verify that files we're managing have correct checksums."""
        
        badChecksums = []
        for filepath in self.filepaths:
            if os.path.exists(filepath):
                checksum = Checksum(filepath, hashtype=self.header.get('hashtype')).get()
                knownChecksum = self.checksums[filepath]
                if checksum != knownChecksum:
                    badChecksums.append([filepath, checksum, knownChecksum])
        return badChecksums

        
    def create(self, hashtype):
        """Create manifest info (checksums, timestamps) for the directory.

        @param hashtype  Type of checksum to use: crc32 or md5.
        """
        
        def callback(checksums, directory, files):
            for file in files:
                path = os.path.join(directory, file)
                if not re.search(".fits(.gz)?$", path):
                    continue
                checksum = Checksum(path, hashtype=hashtype).get()
                manifestLine = path + " " + str(checksum) + "\n"
                checksums[path] = checksum

        self.header = ManifestHeader()
        self.header.set("hashtype", hashtype)
        date = datetime.datetime.now().strftime("%a %Y-%m-%d %H:%M:%S")
        self.header.set("date", date)

        self.checksums = {}
        os.path.walk(self.testdataDir, callback, self.checksums)
        self.filepaths = self.checksums.keys()

        
    def write(self):
        """Write ourself as a manifest file. """
        
        fp = open(self.manifest, 'w')
        fp.write(self.header.write())
        for path in self.filepaths:
            line = path + " " + str(self.checksums[path]) + "\n"
            fp.write(line)
        fp.close()



def verifyManifest(dir, verifyExists=True, verifyChecksum=True, raiseOnFailure=True):
    """Verify the data files in a manifest.

    @param verifyExists     Verify that files exist.
    @param verifyChecksum   Verify the files have the correct checksum
    @param raiseOnFailure   Raise an exception if a file fails.
    """
    
    manifest = Manifest(dir)
    manifest.read()

    msg = ""
    
    missingInputs = []
    if verifyExists:
        missingInputs   = manifest.verifyExists()
    if (len(missingInputs) > 0):
        msg = "Missing input files listed in manifest:\n"
        msg += "\n".join(missingInputs) + "\n"
        
    failedChecksums = []
    if verifyChecksum:
        failedChecksums = manifest.verifyChecksum()
    if len(failedChecksums) > 0:
        msg += "Failed checksums:\n"
        msg += "\n".join(failedChecksums) + "\n"
    if len(msg) > 1 and raiseOnFailure:
        raise Exception(msg)

    return missingInputs, failedChecksums
