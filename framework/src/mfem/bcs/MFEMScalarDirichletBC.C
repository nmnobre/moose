//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMScalarDirichletBC.h"

registerMooseObject("MooseApp", MFEMScalarDirichletBC);

InputParameters
MFEMScalarDirichletBC::validParams()
{
  InputParameters params = MFEMEssentialBC::validParams();
  params.addClassDescription("Applies a Dirichlet condition to a scalar variable.");
  params.addParam<MFEMScalarCoefficientName>(
      "coefficient", "0.", "The coefficient setting the values on the essential boundary");
  return params;
}

MFEMScalarDirichletBC::MFEMScalarDirichletBC(const InputParameters & parameters)
  : MFEMEssentialBC(parameters), _coef(getScalarCoefficient("coefficient"))
{
}

void
MFEMScalarDirichletBC::ApplyBC(mfem::GridFunction & gridfunc)
{
  gridfunc.ProjectBdrCoefficient(_coef, getBoundaryMarkers());
}

void
MFEMScalarDirichletBC::ApplyBC(mfem::ParComplexGridFunction & gridfunc)
{
  if (_real)
  {
    gridfunc.real().SyncMemory(gridfunc);
    gridfunc.real().ProjectBdrCoefficient(_coef, getBoundaryMarkers());
    gridfunc.real().SyncAliasMemory(gridfunc);
  }
  if (_imag)
  {
    gridfunc.imag().SyncMemory(gridfunc);
    gridfunc.imag().ProjectBdrCoefficient(_coef, getBoundaryMarkers());
    gridfunc.imag().SyncAliasMemory(gridfunc);
  }
}

#endif
