# -*- python -*-
#
# Setup our environment
#
import glob, os.path, re
import lsst.SConsUtils as scons

thisPkg    = "testing_pipeQA"
pythonPkg  = "pipeQaLib"
pythonPath = os.path.join("python", "lsst", "testing", "pipeQA")


###############################################################################
# Boilerplate below here
try:
    scons.ConfigureDependentProducts
except AttributeError:
    import lsst.afw.SconsUtils
    scons.ConfigureDependentProducts = lsst.afw.SconsUtils.ConfigureDependentProducts


dependencies = ["pipette"]

#libs = "meas_algorithms ndarray afw daf_base daf_data daf_persistence "
#libs += "pex_logging pex_exceptions pex_policy security boost minuit2 utils wcslib"

env = scons.makeEnv(thisPkg,
                    r"$HeadURL: $",
                    scons.ConfigureDependentProducts(thisPkg))

#env.libs[thisPkg] += env.getlibs(libs) #" ".join(dependencies))

env.thisPkg    = thisPkg
env.pythonPkg  = pythonPkg
env.pythonPath = pythonPath

env.Help("""
Pipeline output testing package
""")


# Build/install things
#SConscript(os.path.join("tests", "SConscript"))
SConscript(os.path.join(pythonPath, "SConscript"))
SConscript(os.path.join("lib", "SConscript"))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"
Alias("install", [])
scons.CleanTree(r"*~ core *.so *.os *.o")

