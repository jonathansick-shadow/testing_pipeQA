#include<cmath>
#include<limits>

#include "lsst/testing/pipeQA/SimRefObject.h"

float const fNaN = std::numeric_limits<float>::quiet_NaN();

void SimRefObject::setMag(double magNew, char filter) {
    
    switch (filter) {
      case 'u':
        _uMag = static_cast<float>(magNew);
        break;
      case 'g':
        _gMag = static_cast<float>(magNew);
        break;
      case 'r':
        _rMag = static_cast<float>(magNew);
        break;
      case 'i':
        _iMag = static_cast<float>(magNew);
        break;
      case 'z':
        _zMag = static_cast<float>(magNew);
        break;
      default:
        break;
    } 
}


void SimRefObject::setFlux(double fluxNew, char filter) {
    if (fluxNew > 0 && ! ::isnan(fluxNew)) {
        setMag(-2.5* ::log(fluxNew)/::log(10.0), filter);
    } else {
        setMag(fNaN, filter);
    }
}


float SimRefObject::getMag(char filter) {
    switch (filter) {
      case 'u':
        return _uMag;
        break;
      case 'g':
        return _gMag;
        break;
      case 'r':
        return _rMag;
        break;
      case 'i':
        return _iMag;
        break;
      case 'z':
        return _zMag;
        break;
      default:
        break;
    } 
    return fNaN;
}

float SimRefObject::getFlux(char filter) {
    return ::pow(10.0, -0.4*getMag(filter));
}
    
template class SimRefObjectSet<SimRefObject>;
