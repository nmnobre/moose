mu = 1.0
epsilon = 1.0
sigma = 20.0
omega = 10.0
kappa_r = 12.7201964951406889525742371916
kappa_i = -7.86151377757423297509831172647

[Mesh]
  type = MFEMMesh
  file = ../mesh/inline-quad.mesh
[]

[Problem]
  type = MFEMProblem
  numeric_type = complex
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
[]

[Variables]
  [u]
    type = MFEMComplexVariable
    fespace = H1FESpace
  []
[]

[Functions]
  [u0_r]
    type = ParsedFunction
    expression = exp(y*${kappa_i})*cos(y*${kappa_r})
  []
  [u0_i]
    type = ParsedFunction
    expression = -exp(y*${kappa_i})*sin(y*${kappa_r})
  []
  [stiffnessCoef]
    type = ParsedFunction
    expression = 1.0/${mu}
  []
  [massCoef]
    type = ParsedFunction
    expression = -${omega}*${omega}*${epsilon}
  []
  [lossCoef]
    type = ParsedFunction
    expression = ${omega}*${sigma}
  []
[]

[BCs]
  [dbc_real]
    type = MFEMScalarDirichletBC
    variable = u
    coefficient = u0_r
    numeric_type = real
  []
  [dbc_imag]
    type = MFEMScalarDirichletBC
    variable = u
    coefficient = u0_i
    numeric_type = imag
  []
[]

[Kernels]
  [diffusion]
    type = MFEMDiffusionKernel
    variable = u
    coefficient = stiffnessCoef
  []
  [mass_real]
    type = MFEMMassKernel
    variable = u
    coefficient = massCoef
    numeric_type = real
  []
  [mass_imag]
    type = MFEMMassKernel
    variable = u
    coefficient = lossCoef
    numeric_type = imag
  []
[]

[Solvers]
  [main]
    type = MFEMSuperLU
  []
[]

[Executioner]
  type = MFEMSteady
  assembly_level = legacy
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/Complex2DQuad
    vtk_format = ASCII
  []
[]
