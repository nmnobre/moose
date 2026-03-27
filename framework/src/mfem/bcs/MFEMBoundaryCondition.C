//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMBoundaryCondition.h"
#include "MFEMProblem.h"
#include "libmesh/ignore_warnings.h"
#include "mfem/miniapps/common/mesh_extras.hpp"
#include "libmesh/restore_warnings.h"

InputParameters
MFEMBoundaryCondition::validParams()
{
  InputParameters params = MFEMObject::validParams();
  params += MFEMBoundaryRestrictable::validParams();

  params.addClassDescription("Base class for applying boundary conditions to MFEM problems.");
  params.registerBase("BoundaryCondition");
  params.registerSystemAttributeName("BoundaryCondition");
  params.addParam<VariableName>("variable", "Variable on which to apply the boundary condition");
  MooseEnum numeric_types("real imag complex", "real");
  params.addParam<MooseEnum>("numeric_type", numeric_types, "Number type used for the bc");
  return params;
}

MFEMBoundaryCondition::MFEMBoundaryCondition(const InputParameters & parameters)
  : MFEMObject(parameters),
    MFEMBoundaryRestrictable(
        parameters, getMFEMProblem().getMFEMVariableMesh(getParam<VariableName>("variable"))),
    _test_var_name(getParam<VariableName>("variable")),
    _real(getParam<MooseEnum>("numeric_type") == "real" ||
          getParam<MooseEnum>("numeric_type") == "complex"),
    _imag(getParam<MooseEnum>("numeric_type") == "imag" ||
          getParam<MooseEnum>("numeric_type") == "complex")
{
  mooseAssert(!_imag || getMFEMProblem().getProblemData().cmplx_gridfunctions.Has(_test_var_name),
              "Requested imag component, but variable is not complex.");
}

#endif
