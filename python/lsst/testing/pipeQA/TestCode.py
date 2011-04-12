# if testing/displayQA is setup, use the www display
# otherwise, just use the default TestSet

# each of these modules must define:
#- TestSet
#- Test
try:
    from lsst.testing.displayQA.TestCode import *
except Exception, e:
    from DefaultTestCode import *

