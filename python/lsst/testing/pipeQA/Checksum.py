import hashlib
import zlib
import re

def hashtypesDefined():
    return ["crc32", "md5"]

class Checksum(object):
    """A class to manage file checksums."""

    def __init__(self, path, bufsize=None, hashtype="crc32"):
        """
        @param path     The file to checksum
        @param bufsize  The buffer to checksum (whole file if set to None)
        @param hashtype The type of checksum to use (crc32 by default)
        """
        
        self.path = path
        self.bufsize = bufsize
        self.hashtype = hashtype

        # open the file and store the requested buffer size worth of data
        fp = open(self.path, 'r')
        if bufsize is not None:
            self.buffer = fp.read(bufsize)
        else:
            self.buffer = fp.read()
        fp.close()

        
    def getMd5(self):
        """Get an md5 checksum for the buffer."""
        md5 = hashlib.md5()
        md5.update(self.buffer)
        return md5.hexdigest()
    
    def getCrc32(self):
        """Get a regular 32 bit checksum for the buffer."""
        crc32 = zlib.crc32(self.buffer)
        return str(crc32)  #md5 returns str

    def get(self):
        """Get the checksum for the hashtype specified at instantiation."""
        if re.search("crc32", self.hashtype):
            return self.getCrc32()
        if re.search("md5", self.hashtype):
            return self.getMd5()
