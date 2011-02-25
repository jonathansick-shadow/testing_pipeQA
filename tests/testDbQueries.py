#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import os
import eups
import lsst.testing.pipeQA as pipeQA

simdir = eups.productDir("obs_lsstSim")
cameraGeomPaf = os.path.join(simdir, "description", "Full_STA_geom.paf")
print cameraGeomPaf

 
foo = pipeQA.ZeropointFigure(cameraGeomPaf, "rplante_DC3b_u_weeklytest_2011_0218_science")
visitId = 85661762
filter  = "r"
foo.fillValues(visitId, filter)
foo.makeFigure("%d %s" % (visitId, filter), doLabel = True)
foo.saveFigure("foo.png")
               
