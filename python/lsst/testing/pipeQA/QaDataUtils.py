import os, re

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
    

