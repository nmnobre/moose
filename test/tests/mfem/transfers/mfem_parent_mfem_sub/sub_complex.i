[Mesh]
  type = MFEMMesh
  file = ../../mesh/square.msh
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

[BCs]
  [bottom]
    type = MFEMScalarDirichletBC
    variable = u
    boundary = 2
    coefficient = 1
    numeric_type = complex
  []
  [top]
    type = MFEMScalarDirichletBC
    variable = u
    boundary = 4
    numeric_type = complex
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = u
    numeric_type = complex
  []
[]

[Solvers]
  [main]
    type = MFEMSuperLU
  []
[]

[Executioner]
  type = MFEMSteady
[]

[MultiApps]
  active = ''
  [subapp]
    type = FullSolveMultiApp
    input_files = parent_complex.i
    execute_on = FINAL
  []
[]

[Transfers]
  active = ''
  [to_sub]
    type = MultiAppMFEMCopyTransfer
    source_variables = u
    variables = u
    to_multi_app = subapp
  []
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/DiffusionSubComplex
    vtk_format = ASCII
  []
[]
