import os, re
import datetime
from lsst.testing.pipeQA.Checksum import Checksum

# could inherit from Dict perhaps
class ManifestHeader(object):

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
        self.entries[key] = value

    def get(self, key):
        return self.entries[key]

    
    def verify(self):

        # check for missing keys
        missingKeys = []
        for k in self.requiredKeys:
            if self.entries[k] is None:
                missingKeys.append(k)
        return missingKeys

    
    def write(self):

        missingKeys = self.verify()
        if missingKeys:
            print "Warning: required manifest header keys missing:\n" + "\n".join(missingKeys)

        outstring = ""
        for k,v in self.entries.items():
            outstring += "# " + k + " " + v + "\n"
        return outstring


    
class Manifest(object):

    def __init__(self, testdataDir):
        self.testdataDir  = testdataDir
        self.checksums    = {}
        self.filepaths    = []
        self.header       = ManifestHeader([])
        
        
    def read(self):
        manifest = os.path.join(self.testdataDir, "manifest")
        missingInputs = []
        failedMd5s = []

        fp = open(manifest, 'r')
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
        return self.header.get(key)

    def verifyExists(self):
        missingInputs = []
        for filepath in self.filepaths:
            if not os.path.exists(filepath):
                missingInputs.append(filepath)
        return missingInputs

    def verifyChecksum(self):
        badChecksums = []
        for filepath in self.filepaths:
            if os.path.exists(filepath):
                checksum = Checksum(filepath, hashtype=self.header.get('hashtype')).get()
                knownChecksum = self.checksums[filepath]
                if checksum != knownChecksum:
                    badChecksums.append([filepath, checksum, knownChecksum])
        return badChecksums

        
    def create(self, hashtype):
        
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
        
        manifest = os.path.join(self.testdataDir, "manifest")
        fp = open(manifest, 'w')
        fp.write(self.header.write())
        for path in self.filepaths:
            line = path + " " + str(self.checksums[path]) + "\n"
            fp.write(line)
        fp.close()


