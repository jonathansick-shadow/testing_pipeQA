#include<vector>



class SimRefObject {
public:
    SimRefObject() :
        _id(0), _isStar(0),
        _ra(0.0), _dec(0.0),
        _uMag(0.0), _gMag(0.0), _rMag(0.0), _iMag(0.0), _zMag(0.0) {}

    SimRefObject(long id, int isStar,
                 double ra, double dec,
                 float uMag, float gMag, float rMag, float iMag, float zMag) :
        _id(id), _isStar(isStar),
        _ra(ra), _dec(dec),
        _uMag(uMag), _gMag(gMag), _rMag(rMag), _iMag(iMag), _zMag(zMag) {}

    long getId() { return _id; }
    void setId(long val) { _id = val; }

    int getIsStar() { return _isStar; }
    void setIsStar(int val) { _isStar = val; }

    double getRa() { return _ra; }
    void setRa(double ra) { _ra = ra; }

    double getDec() { return _dec; }
    void setDec(double dec) { _dec = dec; }

    void setMag(double magNew, char filter);
    void setFlux(double fluxNew, char filter);        
    float getMag(char filter);
    float getFlux(char filter);
    
private:
    long _id;
    int _isStar;
    double _ra, _dec;
    float _uMag, _gMag, _rMag, _iMag, _zMag;
};


template<typename SRO>
class SimRefObjectSet : public std::vector<SRO> {
    SimRefObjectSet() : std::vector<SRO>() {}
    void append(SRO sro) {
        this->push_back(sro);
    }
};

// typedef std::vector<SimRefObject> SimRefObjectSet;
