import hashlib
import zlib
import re

def hashtypesDefined():
    return ["crc32", "md5"]

class Checksum(object):

    def __init__(self, path, bufsize=None, hashtype="crc32"):
        self.path = path
        self.bufsize = bufsize
        self.hashtype = hashtype
        
        fp = open(self.path, 'r')
        if bufsize is not None:
            self.buffer = fp.read(bufsize)
        else:
            self.buffer = fp.read()
        fp.close()

        
    def getMd5(self):
        md5 = hashlib.md5()
        md5.update(self.buffer)
        return md5.hexdigest()
    
    def getCrc32(self):
        crc32 = zlib.crc32(self.buffer)
        return str(crc32)  #md5 returns str

    def get(self):
        if re.search("crc32", self.hashtype):
            return self.getCrc32()
        if re.search("md5", self.hashtype):
            return self.getMd5()
