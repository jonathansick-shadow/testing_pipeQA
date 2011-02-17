# -*- python -*-
#
# Setup our environment
#
import glob, os.path, re
import lsst.SConsUtils as scons

dependencies = ["pipette"]

env = scons.makeEnv("pipetest",
                    r"$HeadURL: $",
                    [
                    ])
env.Help("""
Pipeline output testing package
""")

###############################################################################
# Boilerplate below here

pkg = env["eups_product"]
env.libs[pkg] += env.getlibs(" ".join(dependencies))

# Build/install things
SConscript(os.path.join("tests", "SConscript"))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"
Alias("install", [])
scons.CleanTree(r"*~ core *.so *.os *.o")

