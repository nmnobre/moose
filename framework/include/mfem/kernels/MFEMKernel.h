//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#pragma once

#include "MFEMObject.h"
#include "MFEMContainers.h"
#include "MFEMBlockRestrictable.h"

/**
 * Class to construct an MFEM integrator to apply to the equation system.
 */
class MFEMKernel : public MFEMObject, public MFEMBlockRestrictable
{
public:
  static InputParameters validParams();

  MFEMKernel(const InputParameters & parameters);

  virtual ~MFEMKernel() = default;

  /// Create MFEM integrator to apply to the RHS of the weak form. Ownership managed by the caller.
  virtual mfem::LinearFormIntegrator * createLFIntegrator() { return nullptr; }

  /// Create MFEM integrators to apply to the RHS of the weak form. Ownership managed by the caller.
  virtual std::pair<mfem::LinearFormIntegrator *, mfem::LinearFormIntegrator *>
  createComplexLFIntegrator() final
  {
    return std::make_pair(_real ? createRealLFIntegrator() : nullptr,
                          _imag ? createImagLFIntegrator() : nullptr);
  }

  /// Create MFEM integrator to apply to the LHS of the weak form. Ownership managed by the caller.
  virtual mfem::BilinearFormIntegrator * createBFIntegrator() { return nullptr; }

  /// Create MFEM integrators to apply to the LHS of the weak form. Ownership managed by the caller.
  virtual std::pair<mfem::BilinearFormIntegrator *, mfem::BilinearFormIntegrator *>
  createComplexBFIntegrator() final
  {
    return std::make_pair(_real ? createRealBFIntegrator() : nullptr,
                          _imag ? createImagBFIntegrator() : nullptr);
  }

  /// Create MFEM non-linear integrator to apply to the LHS of the weak form. Ownership managed by the caller.
  virtual mfem::NonlinearFormIntegrator * createNLIntegrator() { return nullptr; }

  /// Get name of the test variable labelling the weak form this kernel is added to
  const VariableName & getTestVariableName() const { return _test_var_name; }

  /// Get name of the trial variable (gridfunction) the kernel acts on.
  /// Defaults to the name of the test variable labelling the weak form.
  virtual const VariableName & getTrialVariableName() const { return _test_var_name; }

protected:
  /// Create real component of the MFEM integrator to apply to the RHS of the weak form.
  virtual mfem::LinearFormIntegrator * createRealLFIntegrator() { return createLFIntegrator(); }

  /// Create imag component of the MFEM integrator to apply to the RHS of the weak form.
  virtual mfem::LinearFormIntegrator * createImagLFIntegrator() { return createLFIntegrator(); }

  /// Create real component of the MFEM integrator to apply to the LHS of the weak form.
  virtual mfem::BilinearFormIntegrator * createRealBFIntegrator() { return createBFIntegrator(); }

  /// Create imag component of the MFEM integrator to apply to the LHS of the weak form.
  virtual mfem::BilinearFormIntegrator * createImagBFIntegrator() { return createBFIntegrator(); }

  /// Name of (the test variable associated with) the weak form that the kernel is applied to.
  const VariableName & _test_var_name;

  /// Whether this is a real kernel contribution
  const bool _real;

  /// Whether this is an imag kernel contribution
  const bool _imag;
};

#endif
