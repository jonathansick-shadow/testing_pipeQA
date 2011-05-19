# if testing/displayQA is setup, use the full QaFigure
# otherwise, just use the default one here (no maps area support)

# each of these modules must define:
#- QaFigure
try:
    from lsst.testing.displayQA.figures import *
except Exception, e:
    print "testing_displayQA not available.  Loading DefaultQaFigure."
    print "Error: ", e
    from DefaultQaFigure import *


