//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "libmesh/ignore_warnings.h"
#include "mfem/miniapps/common/mfem-common.hpp"
#include "libmesh/restore_warnings.h"
#include "MFEMObjectUnitTest.h"
#include "MFEMScalarDirichletBC.h"
#include "MFEMVectorDirichletBC.h"
#include "MFEMVectorNormalDirichletBC.h"
#include "MFEMVectorTangentialDirichletBC.h"

class MFEMEssentialBCTest : public MFEMObjectUnitTest
{
public:
  mfem::common::H1_ParFESpace _scalar_fes;
  mfem::common::H1_ParFESpace _vector_h1_fes;
  mfem::common::ND_ParFESpace _vector_hcurl_fes;
  mfem::common::RT_ParFESpace _vector_hdiv_fes;
  mfem::ParGridFunction _scalar_gridfunc, _vector_h1_gridfunc, _vector_hcurl_gridfunc,
      _vector_hdiv_gridfunc;
  mfem::ParComplexGridFunction _scalar_cgridfunc, _vector_h1_cgridfunc, _vector_hcurl_cgridfunc,
      _vector_hdiv_cgridfunc;

  MFEMEssentialBCTest()
    : MFEMObjectUnitTest("MooseUnitApp"),
      _scalar_fes(_mfem_mesh_ptr->getMFEMParMeshPtr().get(), 1, 3),
      _vector_h1_fes(
          _mfem_mesh_ptr->getMFEMParMeshPtr().get(), 1, 3, mfem::BasisType::GaussLobatto, 3),
      _vector_hcurl_fes(_mfem_mesh_ptr->getMFEMParMeshPtr().get(), 2, 3),
      _vector_hdiv_fes(_mfem_mesh_ptr->getMFEMParMeshPtr().get(), 2, 3),
      _scalar_gridfunc(&_scalar_fes),
      _vector_h1_gridfunc(&_vector_h1_fes),
      _vector_hcurl_gridfunc(&_vector_hcurl_fes),
      _vector_hdiv_gridfunc(&_vector_hdiv_fes),
      _scalar_cgridfunc(&_scalar_fes),
      _vector_h1_cgridfunc(&_vector_h1_fes),
      _vector_hcurl_cgridfunc(&_vector_hcurl_fes),
      _vector_hdiv_cgridfunc(&_vector_hdiv_fes)
  {
    InputParameters func_params = _factory.getValidParams("ParsedFunction");
    func_params.set<std::string>("expression") = "x + y";
    _mfem_problem->addFunction("ParsedFunction", "func1", func_params);
    _mfem_problem->getFunction("func1").initialSetup();
    InputParameters vecfunc_params = _factory.getValidParams("ParsedVectorFunction");
    vecfunc_params.set<std::string>("expression_x") = "x + y";
    vecfunc_params.set<std::string>("expression_y") = "x + y + 1";
    vecfunc_params.set<std::string>("expression_z") = "x + y + 2";
    _mfem_problem->addFunction("ParsedVectorFunction", "func2", vecfunc_params);
    _mfem_problem->getFunction("func2").initialSetup();
    _scalar_gridfunc = 0;
    _vector_h1_gridfunc = 0;
    _vector_hcurl_gridfunc = 0;
    _vector_hdiv_gridfunc = 0;
    _scalar_cgridfunc = std::complex<mfem::real_t>{0, 0};
    _vector_h1_cgridfunc = std::complex<mfem::real_t>{0, 0};
    _vector_hcurl_cgridfunc = std::complex<mfem::real_t>{0, 0};
    _vector_hdiv_cgridfunc = std::complex<mfem::real_t>{0, 0};
    // Register a dummy (Par)GridFunction for the variable the BCs apply to
    auto pgf = std::make_shared<mfem::ParGridFunction>(&_scalar_fes);
    auto pcgf = std::make_shared<mfem::ParComplexGridFunction>(&_scalar_fes);
    _mfem_problem->getProblemData().gridfunctions.Register("test_variable_name", pgf);
    _mfem_problem->getProblemData().cmplx_gridfunctions.Register("test_cmplx_variable_name", pcgf);
  }

  void check_boundary(int /*bound*/,
                      mfem::FiniteElementSpace & fespace,
                      std::function<mfem::real_t(mfem::ElementTransformation *,
                                                 const mfem::IntegrationPoint &)> error_func,
                      mfem::real_t tolerance)
  {
    for (int be = 0; be < _mfem_mesh_ptr->getMFEMParMeshPtr()->GetNBE(); be++)
    {
      mfem::Element * elem = _mfem_mesh_ptr->getMFEMParMeshPtr()->GetBdrElement(be);
      if (elem->GetAttribute() != 1)
        continue;
      mfem::ElementTransformation * trans =
          _mfem_mesh_ptr->getMFEMParMeshPtr()->GetBdrElementTransformation(be);
      const mfem::FiniteElement * fe = fespace.GetBE(be);
      const mfem::IntegrationRule & ir =
          mfem::IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder() + 2);
      mfem::real_t total_error = 0.0;
      for (int j = 0; j < ir.GetNPoints(); j++)
      {
        const mfem::IntegrationPoint point = ir.IntPoint(j);
        trans->SetIntPoint(&point);
        const mfem::real_t error = error_func(trans, point);
        total_error += error * error;
      }
      EXPECT_LT(total_error, tolerance);
    }
  }

  mfem::Vector calc_normal(mfem::ElementTransformation * trans) const
  {
    mfem::Vector normal(3);
    mfem::CalcOrtho(trans->Jacobian(), normal);
    return normal;
  }
};

/**
 * Test MFEMScalarDirichletBC can be constructed from a constant and applied successfully
 */
TEST_F(MFEMEssentialBCTest, MFEMScalarDirichletConstantBC)
{
  // Construct boundary condition
  InputParameters bc_params = _factory.getValidParams("MFEMScalarDirichletBC");
  bc_params.set<VariableName>("variable") = "test_variable_name";
  bc_params.set<MFEMScalarCoefficientName>("coefficient") = "1.";
  bc_params.set<std::vector<BoundaryName>>("boundary") = {"1"};
  auto & essential_bc = addObject<MFEMScalarDirichletBC>("MFEMScalarDirichletBC", "bcr", bc_params);

  EXPECT_EQ(essential_bc.getTrialVariableName(), "test_variable_name");
  EXPECT_EQ(essential_bc.getTestVariableName(), "test_variable_name");

  // Test applying the BC
  essential_bc.ApplyBC(_scalar_gridfunc);

  // Check the correct boundary values have been applied
  mfem::GridFunctionCoefficient scalar_variable(&_scalar_gridfunc);
  check_boundary(
      1,
      _scalar_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      { return std::abs(scalar_variable.Eval(*trans, point) - 1.); },
      1e-8);

  bc_params.set<VariableName>("variable") = "test_cmplx_variable_name";
  bc_params.set<MooseEnum>("numeric_type") = "complex";
  auto & cmplx_essential_bc =
      addObject<MFEMScalarDirichletBC>("MFEMScalarDirichletBC", "bcc", bc_params);

  EXPECT_EQ(cmplx_essential_bc.getTrialVariableName(), "test_cmplx_variable_name");
  EXPECT_EQ(cmplx_essential_bc.getTestVariableName(), "test_cmplx_variable_name");

  // Test applying the BC
  cmplx_essential_bc.ApplyBC(_scalar_cgridfunc);

  // Check the correct boundary values have been applied
  mfem::GridFunctionCoefficient scalar_variable_real(&_scalar_cgridfunc.real());
  mfem::GridFunctionCoefficient scalar_variable_imag(&_scalar_cgridfunc.imag());
  check_boundary(
      1,
      _scalar_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        return std::abs(scalar_variable_real.Eval(*trans, point) - 1.) +
               std::abs(scalar_variable_imag.Eval(*trans, point) - 1.);
      },
      1e-8);
}

/**
 * Test MFEMScalarDirichletBC can be constructed and applied successfully
 */
TEST_F(MFEMEssentialBCTest, MFEMScalarDirichletBC)
{
  // Construct boundary condition
  InputParameters bc_params = _factory.getValidParams("MFEMScalarDirichletBC");
  bc_params.set<VariableName>("variable") = "test_variable_name";
  bc_params.set<MFEMScalarCoefficientName>("coefficient") = "func1";
  bc_params.set<std::vector<BoundaryName>>("boundary") = {"1"};
  auto & essential_bc = addObject<MFEMScalarDirichletBC>("MFEMScalarDirichletBC", "bcr", bc_params);

  EXPECT_EQ(essential_bc.getTrialVariableName(), "test_variable_name");
  EXPECT_EQ(essential_bc.getTestVariableName(), "test_variable_name");

  // Test applying the BC
  essential_bc.ApplyBC(_scalar_gridfunc);

  // Check the correct boundary values have been applied
  mfem::GridFunctionCoefficient scalar_variable(&_scalar_gridfunc);
  mfem::Coefficient & expected(_mfem_problem->getCoefficients().getScalarCoefficient("func1"));
  check_boundary(
      1,
      _scalar_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      { return std::abs(scalar_variable.Eval(*trans, point) - expected.Eval(*trans, point)); },
      1e-8);

  bc_params.set<VariableName>("variable") = "test_cmplx_variable_name";
  bc_params.set<MooseEnum>("numeric_type") = "imag";
  auto & cmplx_essential_bc =
      addObject<MFEMScalarDirichletBC>("MFEMScalarDirichletBC", "bci", bc_params);

  EXPECT_EQ(cmplx_essential_bc.getTrialVariableName(), "test_cmplx_variable_name");
  EXPECT_EQ(cmplx_essential_bc.getTestVariableName(), "test_cmplx_variable_name");

  // Test applying the BC
  cmplx_essential_bc.ApplyBC(_scalar_cgridfunc);

  // Check the correct boundary values have been applied
  mfem::GridFunctionCoefficient scalar_variable_real(&_scalar_cgridfunc.real());
  mfem::GridFunctionCoefficient scalar_variable_imag(&_scalar_cgridfunc.imag());
  check_boundary(
      1,
      _scalar_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        return std::abs(scalar_variable_real.Eval(*trans, point) - 0.) +
               std::abs(scalar_variable_imag.Eval(*trans, point) - expected.Eval(*trans, point));
      },
      1e-8);
}

/**
 * Test MFEMVectorDirichletBC can be constructed from a constant and applied successfully
 */
TEST_F(MFEMEssentialBCTest, MFEMVectorDirichletConstantBC)
{
  // Construct boundary condition
  InputParameters bc_params = _factory.getValidParams("MFEMVectorDirichletBC");
  bc_params.set<VariableName>("variable") = "test_variable_name";
  bc_params.set<MFEMVectorCoefficientName>("vector_coefficient") = "1. 2. 3.";
  bc_params.set<std::vector<BoundaryName>>("boundary") = {"1"};
  auto & essential_bc = addObject<MFEMVectorDirichletBC>("MFEMVectorDirichletBC", "bc1", bc_params);

  EXPECT_EQ(essential_bc.getTrialVariableName(), "test_variable_name");
  EXPECT_EQ(essential_bc.getTestVariableName(), "test_variable_name");

  // Test applying the BC
  essential_bc.ApplyBC(_vector_h1_gridfunc);

  // Check the correct boundary values have been applied
  mfem::VectorGridFunctionCoefficient variable(&_vector_h1_gridfunc);
  mfem::Vector expected({1., 2., 3.});
  check_boundary(
      1,
      _vector_h1_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        mfem::Vector actual(3);
        variable.Eval(actual, *trans, point);
        actual -= expected;
        return actual.Norml2();
      },
      1e-8);
}

/**
 * Test MFEMVectorDirichletBC can be constructed and applied successfully
 */
TEST_F(MFEMEssentialBCTest, MFEMVectorDirichletBC)
{
  // Construct boundary condition
  InputParameters bc_params = _factory.getValidParams("MFEMVectorDirichletBC");
  bc_params.set<VariableName>("variable") = "test_variable_name";
  bc_params.set<MFEMVectorCoefficientName>("vector_coefficient") = "func2";
  bc_params.set<std::vector<BoundaryName>>("boundary") = {"1"};
  auto & essential_bc = addObject<MFEMVectorDirichletBC>("MFEMVectorDirichletBC", "bc1", bc_params);

  EXPECT_EQ(essential_bc.getTrialVariableName(), "test_variable_name");
  EXPECT_EQ(essential_bc.getTestVariableName(), "test_variable_name");

  // Test applying the BC
  essential_bc.ApplyBC(_vector_h1_gridfunc);

  // Check the correct boundary values have been applied
  mfem::VectorGridFunctionCoefficient variable(&_vector_h1_gridfunc);
  mfem::VectorCoefficient & function(
      _mfem_problem->getCoefficients().getVectorCoefficient("func2"));
  check_boundary(
      1,
      _vector_h1_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        mfem::Vector actual(3), expected(3);
        variable.Eval(actual, *trans, point);
        function.Eval(expected, *trans, point);
        actual -= expected;
        return actual.Norml2();
      },
      1e-8);
}

/**
 * Test MFEMVectorNormalDirichletBC can be constructed from a constant and applied
 * successfully
 */
TEST_F(MFEMEssentialBCTest, MFEMVectorNormalDirichletConstantBC)
{
  // Construct boundary condition
  InputParameters bc_params = _factory.getValidParams("MFEMVectorNormalDirichletBC");
  bc_params.set<VariableName>("variable") = "test_variable_name";
  bc_params.set<MFEMVectorCoefficientName>("vector_coefficient") = "1. 2. 3.";
  bc_params.set<std::vector<BoundaryName>>("boundary") = {"1"};
  auto & essential_bc =
      addObject<MFEMVectorNormalDirichletBC>("MFEMVectorNormalDirichletBC", "bcr", bc_params);

  EXPECT_EQ(essential_bc.getTrialVariableName(), "test_variable_name");
  EXPECT_EQ(essential_bc.getTestVariableName(), "test_variable_name");

  // Test applying the BC
  essential_bc.ApplyBC(_vector_hdiv_gridfunc);

  // Check the correct boundary values have been applied
  mfem::VectorGridFunctionCoefficient variable(&_vector_hdiv_gridfunc);
  mfem::Vector expected({1., 2., 3.});
  check_boundary(
      1,
      _vector_hdiv_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        mfem::Vector actual(3), normal = calc_normal(trans);
        variable.Eval(actual, *trans, point);
        actual -= expected;
        // (actual - expected) should be perpendicular to the normal.
        return std::abs(normal * actual);
      },
      1e-8);

  bc_params.set<VariableName>("variable") = "test_cmplx_variable_name";
  bc_params.set<MooseEnum>("numeric_type") = "complex";
  auto & cmplx_essential_bc =
      addObject<MFEMVectorNormalDirichletBC>("MFEMVectorNormalDirichletBC", "bcc", bc_params);

  EXPECT_EQ(cmplx_essential_bc.getTrialVariableName(), "test_cmplx_variable_name");
  EXPECT_EQ(cmplx_essential_bc.getTestVariableName(), "test_cmplx_variable_name");

  // Test applying the BC
  cmplx_essential_bc.ApplyBC(_vector_hdiv_cgridfunc);

  // Check the correct boundary values have been applied
  mfem::VectorGridFunctionCoefficient variable_real(&_vector_hdiv_cgridfunc.real());
  mfem::VectorGridFunctionCoefficient variable_imag(&_vector_hdiv_cgridfunc.imag());
  check_boundary(
      1,
      _vector_hdiv_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        mfem::Vector actual_real(3), actual_imag(3), normal = calc_normal(trans);
        variable_real.Eval(actual_real, *trans, point);
        variable_imag.Eval(actual_imag, *trans, point);
        actual_real -= expected;
        actual_imag -= expected;
        // (actual - expected) should be perpendicular to the normal.
        return std::abs(normal * actual_real) + std::abs(normal * actual_imag);
      },
      1e-8);
}

/**
 * Test MFEMVectorNormalDirichletBC can be constructed and applied successfully
 */
TEST_F(MFEMEssentialBCTest, MFEMVectorNormalDirichletBC)
{
  // Construct boundary condition
  InputParameters bc_params = _factory.getValidParams("MFEMVectorNormalDirichletBC");
  bc_params.set<VariableName>("variable") = "test_variable_name";
  bc_params.set<MFEMVectorCoefficientName>("vector_coefficient") = "func2";
  bc_params.set<std::vector<BoundaryName>>("boundary") = {"1"};
  auto & essential_bc =
      addObject<MFEMVectorNormalDirichletBC>("MFEMVectorNormalDirichletBC", "bcr", bc_params);

  EXPECT_EQ(essential_bc.getTrialVariableName(), "test_variable_name");
  EXPECT_EQ(essential_bc.getTestVariableName(), "test_variable_name");

  // Test applying the BC
  essential_bc.ApplyBC(_vector_hdiv_gridfunc);

  // Check the correct boundary values have been applied
  mfem::VectorGridFunctionCoefficient variable(&_vector_hdiv_gridfunc);
  mfem::VectorCoefficient & function(
      _mfem_problem->getCoefficients().getVectorCoefficient("func2"));
  check_boundary(
      1,
      _vector_hdiv_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        mfem::Vector actual(3), expected(3), normal = calc_normal(trans);
        variable.Eval(actual, *trans, point);
        function.Eval(expected, *trans, point);
        actual -= expected;
        // (actual - expected) should be perpendicular to the normal.
        return std::abs(normal * actual);
      },
      1e-8);

  bc_params.set<VariableName>("variable") = "test_cmplx_variable_name";
  bc_params.set<MooseEnum>("numeric_type") = "complex";
  auto & cmplx_essential_bc =
      addObject<MFEMVectorNormalDirichletBC>("MFEMVectorNormalDirichletBC", "bcc", bc_params);

  EXPECT_EQ(cmplx_essential_bc.getTrialVariableName(), "test_cmplx_variable_name");
  EXPECT_EQ(cmplx_essential_bc.getTestVariableName(), "test_cmplx_variable_name");

  // Test applying the BC
  cmplx_essential_bc.ApplyBC(_vector_hdiv_cgridfunc);

  // Check the correct boundary values have been applied
  mfem::VectorGridFunctionCoefficient variable_real(&_vector_hdiv_cgridfunc.real());
  mfem::VectorGridFunctionCoefficient variable_imag(&_vector_hdiv_cgridfunc.imag());
  check_boundary(
      1,
      _vector_hdiv_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        mfem::Vector actual_real(3), actual_imag(3), expected(3), normal = calc_normal(trans);
        variable_real.Eval(actual_real, *trans, point);
        variable_imag.Eval(actual_imag, *trans, point);
        function.Eval(expected, *trans, point);
        actual_real -= expected;
        actual_imag -= expected;
        // (actual - expected) should be perpendicular to the normal.
        return std::abs(normal * actual_real) + std::abs(normal * actual_imag);
      },
      1e-8);
}

/**
 * Test MFEMVectorTangentialDirichletBC can be constructed from a constant and applied
 * successfully
 */
TEST_F(MFEMEssentialBCTest, MFEMVectorTangentialDirichletConstantBC)
{
  // Construct boundary condition
  InputParameters bc_params = _factory.getValidParams("MFEMVectorTangentialDirichletBC");
  bc_params.set<VariableName>("variable") = "test_variable_name";
  bc_params.set<MFEMVectorCoefficientName>("vector_coefficient") = "1. 2. 3.";
  bc_params.set<std::vector<BoundaryName>>("boundary") = {"1"};
  auto & essential_bc = addObject<MFEMVectorTangentialDirichletBC>(
      "MFEMVectorTangentialDirichletBC", "bcr", bc_params);

  EXPECT_EQ(essential_bc.getTrialVariableName(), "test_variable_name");
  EXPECT_EQ(essential_bc.getTestVariableName(), "test_variable_name");

  // Test applying the BC
  essential_bc.ApplyBC(_vector_hcurl_gridfunc);
  // Check the correct boundary values have been applied
  mfem::VectorGridFunctionCoefficient variable(&_vector_hcurl_gridfunc);
  mfem::Vector expected({1., 2., 3.});
  check_boundary(
      1,
      _vector_hcurl_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        mfem::Vector actual(3), normal = calc_normal(trans), cross_prod(3);
        variable.Eval(actual, *trans, point);
        actual -= expected;
        // (actual - expected) should be parallel to the normal.
        normal.cross3D(actual, cross_prod);
        return cross_prod.Norml2();
      },
      1e-8);

  bc_params.set<VariableName>("variable") = "test_cmplx_variable_name";
  bc_params.set<MooseEnum>("numeric_type") = "complex";
  auto & cmplx_essential_bc = addObject<MFEMVectorTangentialDirichletBC>(
      "MFEMVectorTangentialDirichletBC", "bcc", bc_params);

  EXPECT_EQ(cmplx_essential_bc.getTrialVariableName(), "test_cmplx_variable_name");
  EXPECT_EQ(cmplx_essential_bc.getTestVariableName(), "test_cmplx_variable_name");

  // Test applying the BC
  cmplx_essential_bc.ApplyBC(_vector_hcurl_cgridfunc);

  // Check the correct boundary values have been applied
  mfem::VectorGridFunctionCoefficient variable_real(&_vector_hcurl_cgridfunc.real());
  mfem::VectorGridFunctionCoefficient variable_imag(&_vector_hcurl_cgridfunc.imag());
  check_boundary(
      1,
      _vector_hcurl_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        mfem::Vector actual_real(3), actual_imag(3), normal = calc_normal(trans),
                                                     cross_prod_real(3), cross_prod_imag(3);
        variable_real.Eval(actual_real, *trans, point);
        variable_imag.Eval(actual_imag, *trans, point);
        actual_real -= expected;
        actual_imag -= expected;
        // (actual - expected) should be parallel to the normal.
        normal.cross3D(actual_real, cross_prod_real);
        normal.cross3D(actual_imag, cross_prod_imag);
        return cross_prod_real.Norml2() + cross_prod_imag.Norml2();
      },
      1e-8);
}

/**
 * Test MFEMVectorTangentialDirichletBC can be constructed and applied successfully
 */
TEST_F(MFEMEssentialBCTest, MFEMVectorTangentialDirichletBC)
{
  // Construct boundary condition
  InputParameters bc_params = _factory.getValidParams("MFEMVectorTangentialDirichletBC");
  bc_params.set<VariableName>("variable") = "test_variable_name";
  bc_params.set<MFEMVectorCoefficientName>("vector_coefficient") = "func2";
  bc_params.set<std::vector<BoundaryName>>("boundary") = {"1"};
  auto & essential_bc = addObject<MFEMVectorTangentialDirichletBC>(
      "MFEMVectorTangentialDirichletBC", "bcr", bc_params);

  EXPECT_EQ(essential_bc.getTrialVariableName(), "test_variable_name");
  EXPECT_EQ(essential_bc.getTestVariableName(), "test_variable_name");

  // Test applying the BC
  essential_bc.ApplyBC(_vector_hcurl_gridfunc);

  // Check the correct boundary values have been applied
  mfem::VectorGridFunctionCoefficient variable(&_vector_hcurl_gridfunc);
  mfem::VectorCoefficient & function(
      _mfem_problem->getCoefficients().getVectorCoefficient("func2"));
  check_boundary(
      1,
      _vector_hcurl_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        mfem::Vector actual(3), expected(3), normal = calc_normal(trans), cross_prod(3);
        variable.Eval(actual, *trans, point);
        function.Eval(expected, *trans, point);
        actual -= expected;
        // (actual - expected) should be parallel to the normal.
        normal.cross3D(actual, cross_prod);
        return cross_prod.Norml2();
      },
      1e-8);

  bc_params.set<VariableName>("variable") = "test_cmplx_variable_name";
  bc_params.set<MooseEnum>("numeric_type") = "complex";
  auto & cmplx_essential_bc = addObject<MFEMVectorTangentialDirichletBC>(
      "MFEMVectorTangentialDirichletBC", "bcc", bc_params);

  EXPECT_EQ(cmplx_essential_bc.getTrialVariableName(), "test_cmplx_variable_name");
  EXPECT_EQ(cmplx_essential_bc.getTestVariableName(), "test_cmplx_variable_name");

  // Test applying the BC
  cmplx_essential_bc.ApplyBC(_vector_hcurl_cgridfunc);

  // Check the correct boundary values have been applied
  mfem::VectorGridFunctionCoefficient variable_real(&_vector_hcurl_cgridfunc.real());
  mfem::VectorGridFunctionCoefficient variable_imag(&_vector_hcurl_cgridfunc.imag());
  check_boundary(
      1,
      _vector_hcurl_fes,
      [&](mfem::ElementTransformation * trans, const mfem::IntegrationPoint & point)
      {
        mfem::Vector actual_real(3), actual_imag(3), expected(3),
            normal = calc_normal(trans), cross_prod_real(3), cross_prod_imag(3);
        variable_real.Eval(actual_real, *trans, point);
        variable_imag.Eval(actual_imag, *trans, point);
        function.Eval(expected, *trans, point);
        actual_real -= expected;
        actual_imag -= expected;
        // (actual - expected) should be parallel to the normal.
        normal.cross3D(actual_real, cross_prod_real);
        normal.cross3D(actual_imag, cross_prod_imag);
        return cross_prod_real.Norml2() + cross_prod_imag.Norml2();
      },
      1e-8);
}

#endif
