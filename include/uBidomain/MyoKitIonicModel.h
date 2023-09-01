#ifndef MYOKITIONICMODEL_H
#define MYOKITIONICMODEL_H

#include "IonicModel.h"
#include "Parameters.h"

namespace uBidomain{
namespace MyoKit{
    
class MyoKitIonicModel: public IonicModel
{
public:
    MyoKitIonicModel();
    MyoKitIonicModel(const Parameters& param_);
    virtual ~MyoKitIonicModel(){};
    
    virtual Vector f(const Real t, const Vector& y, bool Iapp_yes=true) const;        
    void get_y0(Vector& y0) const;
    
private:
    void init();
};

}
}
#endif /* MYOKITIONICMODEL_H */

