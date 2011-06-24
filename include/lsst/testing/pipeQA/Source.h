#include<vector>

class Source {
public:
    Source() {}

    long getId() { return _id; }
    void setId(long val) { _id = val; }

    int getFlagForDetection() { return _flagForDetection; }
    void setFlagForDetection(int val) { _flagForDetection = val; }

    double getRa() { return _ra; }
    void setRa(double ra) { _ra = ra; }

    double getDec() { return _dec; }
    void setDec(double dec) { _dec = dec; }

    void setPhotometry(void *) {}
    void setAstrometry(void *) {}
    void setShape(void*) {}
        
    float getXAstrom() { return _xAstrom; }
    void setXAstrom(float val) { _xAstrom = val; }
    float getYAstrom() { return _yAstrom; }
    void setYAstrom(float val) { _yAstrom = val; }
    
    float getPsfFlux() { return _psfFlux; }
    void setPsfFlux(float val) { _psfFlux = val; }
    float getPsfFluxErr() { return _psfFluxErr; }
    void setPsfFluxErr(float val) { _psfFluxErr = val; }

    float getApFlux() { return _apFlux; }
    void setApFlux(float val) { _apFlux = val; }
    float getApFluxErr() { return _apFluxErr; }
    void setApFluxErr(float val) { _apFluxErr = val; }
    
    float getModelFlux() { return _modelFlux; }
    void setModelFlux(float val) { _modelFlux = val; }
    float getModelFluxErr() { return 0.0; }
    void setModelFluxErr(float val) {}

    float getInstFlux() { return _instFlux; }
    void setInstFlux(float val) { _instFlux = val; }
    float getInstFluxErr() { return 0.0; }
    void setInstFluxErr(float val) { }
    
    float getIxx() { return _ixx; }
    void setIxx(float val) { _ixx = val; }
    float getIxxErr() { return 0.0; }
    void setIxxErr(float val) {}
    float getIyy() { return _iyy; }
    void setIyy(float val) { _iyy = val; }
    float getIyyErr() { return 0.0; }
    void setIyyErr(float val) {}
    float getIxy() { return _ixy; }
    void setIxy(float val) { _ixy = val; }
    float getIxyErr() { return 0.0; }
    void setIxyErr(float val) {}
    

    float getPsfIxx() { return 0.0; }
    void setPsfIxx(float val) {}
    float getPsfIxxErr() { return 0.0; }
    void setPsfIxxErr(float val) {}
    float getPsfIyy() { return 0.0; }
    void setPsfIyy(float val) {}
    float getPsfIyyErr() { return 0.0; }
    void setPsfIyyErr(float val) {}
    float getPsfIxy() { return 0.0; }
    void setPsfIxy(float val) {}
    float getPsfIxyErr() { return 0.0; }
    void setPsfIxyErr(float val) {}
    
    float getResolution() { return 0.0; }
    void setResolution(float val) {}
    
    float getE1() { return 0.0   ; }
    void setE1(float val) {}
    float getE1Err() { return 0.0; }
    void setE1Err(float val) {}
    float getE2() { return 0.0   ; }
    void setE2(float val) {}
    float getE2Err() { return 0.0; }
    void setE2Err(float val) {}
    
    float getShear1() { return 0.0; }
    void setShear1(float val) {}
    float getShear1Err() { return 0.0; }
    void setShear1Err(float val) {}
    float getShear2() { return 0.0; }
    void setShear2(float val) {}
    float getShear2Err() { return 0.0; }
    void setShear2Err(float val) {}
    
    float getSigma() { return 0.0; }
    void setSigma(float val) {}
    float getSigmaErr() { return 0.0; }
    void setSigmaErr(float val) {}
    
private:
    long _id;
    int _flagForDetection;
    double _ra, _dec;
    float _xAstrom, _yAstrom;
    float _psfFlux, _psfFluxErr;
    float _apFlux, _apFluxErr;
    float _instFlux, _modelFlux;
    float _ixx, _iyy, _ixy;
};



template<typename SOURCE>
class SourceSet : public std::vector<SOURCE> {
    SourceSet() : std::vector<SOURCE>() {}
    void append(SOURCE src) {
        this->push_back(src);
    }
};
