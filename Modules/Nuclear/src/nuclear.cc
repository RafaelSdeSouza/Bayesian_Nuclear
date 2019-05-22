/*
 *  Copyright (C) 2017 Rafael S. de Souza <drsouza@ad.unc.edu>
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License v3.0.
 *
 *
 */
#include <module/Module.h>
#include <functions/sfactorTdn.h>
#include <functions/sfactor3Hedp.h>
#include <functions/sfactor3Hedpx.h>
#include <functions/sfactorTdn2.h>
#include <functions/sfactor3Hdnx.h>
#include <functions/sigma7Benpx.h>
#include <functions/sfactor7Lipnx.h>

using std::vector;

namespace jags {
namespace nuclear {

class NUCLEARModule : public Module {
  public:
    NUCLEARModule();
    ~NUCLEARModule();
};

NUCLEARModule::NUCLEARModule() : Module("nuclear"){
  //load functions
     insert(new sfactorTdn);
     insert(new sfactor3Hedp);
     insert(new sfactor3Hedpx);
     insert(new sfactorTdn2);
     insert(new sfactor3Hdnx);
     insert(new  sigma7Benpx);
     insert(new sfactor7Lipnx);
}


NUCLEARModule::~NUCLEARModule() 
{
  vector<Function*> const &fvec = functions();
  for (unsigned int i = 0; i < fvec.size(); ++i) {
    delete fvec[i];
  }
//  vector<Distribution*> const &dvec = distributions();
//  for (unsigned int i = 0; i < dvec.size(); ++i) {
//    delete dvec[i];
//  }
}

} // namespace nuclear
} // namespace jags

jags::nuclear::NUCLEARModule _nuclear_module;
