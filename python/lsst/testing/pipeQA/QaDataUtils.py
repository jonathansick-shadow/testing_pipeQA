import os, re
import math

import lsst.afw.image   as afwImage
import lsst.afw.coord   as afwCoord
import lsst.pex.logging  as pexLog
import lsst.pex.policy  as pexPolicy
import lsst.meas.astrom as measAstrom
import lsst.meas.algorithms.utils as maUtils

def findDataInTestbed(label, raiseOnFailure=True):
    """Scan TESTBED_PATH directories to find a testbed dataset with a given name"""
    
    # If label==None -> use TESTBOT_DIR (_DIR, not _PATH ... only one dataset can be run)
    if re.search("^testBot", label):
        testbotDir = os.getenv("TESTBOT_DIR")
        testbedDir, testdataDir  = os.path.split(testbotDir)
        return testbedDir, testbotDir

    
    # otherwise, get a specific test-data set from one of the testbed directories
    testbedPath = os.getenv("TESTBED_PATH")
    if testbedPath is not None:
        testbedDirs = testbedPath.split(":")
    else:
        raise Exception("Must specify environment variable TESTBED_PATH.")
    
    #############################
    # find the label in the testbed path
    testbedDir = None
    testdataDir = None
    for tbDir in testbedDirs:
        path = os.path.join(tbDir, label)
        if os.path.exists(path):
            testbedDir = tbDir
            testdataDir = path
                
    if raiseOnFailure and (testbedDir is None):
        msg = "Testbed %s was not found in any of the TESTBED_PATH directories:\n" % (label)
        msg += "\n".join(testbedDirs) + "\n"
        raise Exception(msg)
    
    return testbedDir, testdataDir



def getSourceSetNameList():
    
    accessors = [
	["Id",                           "objectId"                      ],
	["Ra",                           "ra",                           ],
	["Dec",                          "decl",                         ],
	["RaErrForWcs",                  "raSigmaForWcs",                ],
	["DecErrForWcs",                 "declSigmaForWcs",              ],
	["RaErrForDetection",            "raSigmaForDetection",          ],
	["DecErrForDetection",           "declSigmaForDetection",        ],
	["XFlux",                        "xFlux",                        ],
	["XFluxErr",                     "xFluxSigma",                   ],
	["YFlux",                        "yFlux",                        ],
	["YFluxErr",                     "yFluxSigma",                   ],
	["RaFlux",                       "raFlux",                       ],
	["RaFluxErr",                    "raFluxSigma",                  ],
	["DecFlux",                      "declFlux",                     ],
	["DecFluxErr",                   "declFluxSigma",                ],
	["XPeak",                        "xPeak",                        ],
	["YPeak",                        "yPeak",                        ],
	["RaPeak",                       "raPeak",                       ],
	["DecPeak",                      "declPeak",                     ],
	["XAstrom",                      "xAstrom",                      ],
	["XAstromErr",                   "xAstromSigma",                 ],
	["YAstrom",                      "yAstrom",                      ],
	["YAstromErr",                   "yAstromSigma",                 ],
	["RaAstrom",                     "raAstrom",                     ],
	["RaAstromErr",                  "raAstromSigma",                ],
	["DecAstrom",                    "declAstrom",                   ],
	["DecAstromErr",                 "declAstromSigma",              ],
	["TaiMidPoint",                  "taiMidPoint",                  ],
	["TaiRange",                     "taiRange",                     ],
	["PsfFlux",                      "psfFlux",                      ],
	["PsfFluxErr",                   "psfFluxSigma",                 ],
	["ApFlux",                       "apFlux",                       ],
	["ApFluxErr",                    "apFluxSigma",                  ],
	["ModelFlux",                    "modelFlux",                    ],
	["ModelFluxErr",                 "modelFluxSigma",               ],
	["InstFlux",                     "instFlux",                     ],
	["InstFluxErr",                  "instFluxSigma",                ],
	["NonGrayCorrFlux",              "nonGrayCorrFlux",              ],
	["NonGrayCorrFluxErr",           "nonGrayCorrFluxSigma",         ],
	["AtmCorrFlux",                  "atmCorrFlux",                  ],
	["AtmCorrFluxErr",               "atmCorrFluxSigma",             ],
	["ApDia",                        "apDia",                        ],
	["Ixx",                          "ixx",                          ],
	["IxxErr",                       "ixxSigma",                     ],
	["Iyy",                          "iyy",                          ],
	["IyyErr",                       "iyySigma",                     ],
	["Ixy",                          "ixy",                          ],
	["IxyErr",                       "ixySigma",                     ],
	["PsfIxx",                       "psfIxx",                       ],
	["PsfIxxErr",                    "psfIxxSigma",                  ],
	["PsfIyy",                       "psfIyy",                       ],
	["PsfIyyErr",                    "psfIyySigma",                  ],
	["PsfIxy",                       "psfIxy",                       ],
	["PsfIxyErr",                    "psfIxySigma",                  ],
	["Resolution",                   "resolution_SG",                ],
	["E1",                           "e1_SG",                        ],
	["E1Err",                        "e1_SG_Sigma",                  ],
	["E2",                           "e2_SG",                        ],
	["E2Err",                        "e2_SG_Sigma",                  ],
	["Shear1",                       "shear1_SG",                    ],
	["Shear1Err",                    "shear1_SG_Sigma",              ],
	["Shear2",                       "shear2_SG",                    ],
	["Shear2Err",                    "shear2_SG_Sigma",              ],
	["Sigma",                        "sourceWidth_SG",               ],
	["SigmaErr",                     "sourceWidth_SG_Sigma",         ],
	#["ShapeStatus",                  "shapeStatus",                  ],
	["Snr",                          "snr",                          ],
	["Chi2",                         "chi2",                         ],
	["FlagForAssociation",           "flagForAssociation",           ],
	["FlagForDetection",             "flagForDetection",             ],
	["FlagForWcs",                   "flagForWcs",                   ],
	]
    return accessors

def getSourceSetAccessors():
    return zip(*getSourceSetNameList())[0]
def getSourceSetDbNames():
    return zip(*getSourceSetNameList())[1]
    

def getSceNameList():
    nameList = [
	['scienceCcdExposureId', 'scienceCcdExposureId' ],
	['visit',                'visit'                ],
	['raft',                 'raft'                 ],
	['raftName',             'raftName'             ],
	['ccd',                  'ccd'                  ],
	['ccdName',              'ccdName'              ],
	['filterId',             'filterId'             ],
	['filterName',           'filterName'           ],
	['ra',                   'ra'                   ],
	['decl',                 'decl'                 ],
	['equinox',              'equinox'              ],
	['raDeSys',              'raDeSys'              ],
	['ctype1',               'ctype1'               ],
	['ctype2',               'ctype2'               ],
	['crpix1',               'crpix1'               ],
	['crpix2',               'crpix2'               ],
	['crval1',               'crval1'               ],
	['crval2',               'crval2'               ],
	['cd1_1',                'cd1_1'                ],
	['cd1_2',                'cd1_2'                ],
	['cd2_1',                'cd2_1'                ],
	['cd2_2',                'cd2_2'                ],
	['llcRa',                'llcRa'                ],
	['llcDecl',              'llcDecl'              ],
	['ulcRa',                'ulcRa'                ],
	['ulcDecl',              'ulcDecl'              ],
	['urcRa',                'urcRa'                ],
	['urcDecl',              'urcDecl'              ],
	['lrcRa',                'lrcRa'                ],
	['lrcDecl',              'lrcDecl'              ],
	['taiMjd',               'taiMjd'               ],
	['obsStart',             'obsStart'             ],
	['expMidpt',             'expMidpt'             ],
	['expTime',              'expTime'              ],
	['nCombine',             'nCombine'             ],
	['binX',                 'binX'                 ],
	['binY',                 'binY'                 ],
	['readNoise',            'readNoise'            ],
	['saturationLimit',      'saturationLimit'      ],
	['gainEff',              'gainEff'              ],
	['fluxMag0',             'fluxMag0'             ],
	['fluxMag0Sigma',        'fluxMag0Sigma'        ],
	['fwhm',                 'fwhm'                 ],
	]
    return nameList

def getSceDbNames():
    return zip(*getSceNameList())[1]




def getCalibObjects(butler, filterName, dataId):
    """
    A version of getCalibObjectsImpl that isn't a class method, for use by other code
    
    @param useOutputSrc             # use fluxes from the "src", not "icSrc"
    """

    log = pexLog.Log.getDefaultLog()
    
    psources = butler.get('icSrc', dataId)
    pmatches = butler.get('icMatch', dataId)
    calexp_md = butler.get('calexp_md', dataId)

    wcs = afwImage.makeWcs(calexp_md)
    W, H = calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2")
    calib = afwImage.Calib(calexp_md)

    matches = pmatches.getSourceMatches()
    sources = psources.getSources()

    anid = pmatches.getSourceMatchMetadata().getInt('ANINDID')

    del psources; del pmatches; del calexp_md # cleanup

    useOutputSrc = False
    if useOutputSrc:
        srcs = butler.get('src', dataId).getSources()
        import lsst.afw.detection as afwDetect
        pmMatch = afwDetect.matchXy(sources, srcs, 1.0, True)
        for icSrc, src, d in pmMatch:
            icSrc.setPsfFlux(src.getPsfFlux())

    # ref sources
    xc, yc = 0.5*W, 0.5*H
    radec = wcs.pixelToSky(xc, yc)
    ra = radec.getLongitude(afwCoord.DEGREES)
    dec = radec.getLatitude(afwCoord.DEGREES)
    radius = wcs.pixelScale()*math.hypot(xc, yc)*1.1

    pol = pexPolicy.Policy()
    pol.set('matchThreshold', 30)
    solver = measAstrom.createSolver(pol, log)
    idName = 'id'

    X = solver.getCatalogue(ra, dec, radius, filterName, idName, anid)
    refsources = X.refsources
    inds = X.inds

    referrs, stargal = None, None
    cols = solver.getTagAlongColumns(anid)
    colnames = [c.name for c in cols]

    col = 'starnotgal'
    if col in colnames:
        stargal1 = solver.getTagAlongBool(anid, col, inds)
        stargal = []
        for i in range(len(stargal1)):
            stargal.append(stargal1[i])

    fdict = maUtils.getDetectionFlags()

    keepref = []
    keepi = []
    for i in xrange(len(refsources)):
        x, y = wcs.skyToPixel(refsources[i].getRa(), refsources[i].getDec())
        if x < 0 or y < 0 or x > W or y > H:
            continue
        refsources[i].setXAstrom(x)
        refsources[i].setYAstrom(y)
        if stargal[i]:
            refsources[i].setFlagForDetection(refsources[i].getFlagForDetection() | fdict["STAR"])
        keepref.append(refsources[i])
        keepi.append(i)

    refsources = keepref

    if referrs is not None:
        referrs = [referrs[i] for i in keepi]
    if stargal is not None:
        stargal = [stargal[i] for i in keepi]

    measAstrom.joinMatchList(matches, refsources, first=True, log=log)
    measAstrom.joinMatchList(matches, sources, first=False, log=log)

    return matches, calib, refsources
