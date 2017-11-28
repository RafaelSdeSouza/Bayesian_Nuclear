/*
 *  Copyright (C) 2017 Rafael S. de Souza <drsouza@ad.unc.edu>
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License v3.0.
 *
 *
 */
#include <module/Module.h>
#include <functions/sfactor3Hedp.h>
#include <functions/sfactorTdn.h>
#include <functions/sfactorTdn_5p.h>

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
  insert(new sfactor3Hedp);
  insert(new sfactorTdn);
  insert(new sfactorTdn_5p);
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
