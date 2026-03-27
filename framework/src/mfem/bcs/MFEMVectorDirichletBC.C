//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMVectorDirichletBC.h"

registerMooseObject("MooseApp", MFEMVectorDirichletBC);

InputParameters
MFEMVectorDirichletBC::validParams()
{
  InputParameters params = MFEMVectorDirichletBCBase::validParams();
  params.addClassDescription(
      "Applies a Dirichlet condition to all components of a vector variable.");
  return params;
}

MFEMVectorDirichletBC::MFEMVectorDirichletBC(const InputParameters & parameters)
  : MFEMVectorDirichletBCBase(parameters)
{
}

void
MFEMVectorDirichletBC::ApplyBC(mfem::GridFunction & gridfunc)
{
  gridfunc.ProjectBdrCoefficient(_vec_coef, getBoundaryMarkers());
}

void
MFEMVectorDirichletBC::ApplyBC(mfem::ParComplexGridFunction & gridfunc)
{
  if (_real)
  {
    gridfunc.real().SyncMemory(gridfunc);
    gridfunc.real().ProjectBdrCoefficient(_vec_coef, getBoundaryMarkers());
    gridfunc.real().SyncAliasMemory(gridfunc);
  }
  if (_imag)
  {
    gridfunc.imag().SyncMemory(gridfunc);
    gridfunc.imag().ProjectBdrCoefficient(_vec_coef, getBoundaryMarkers());
    gridfunc.imag().SyncAliasMemory(gridfunc);
  }
}

#endif
