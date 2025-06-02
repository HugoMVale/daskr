project(
  'daskr',
  ['fortran', 'c'],
  version: '0.0.1',
  meson_version: '>=1.1',
  license: 'MIT',
  default_options: [
  'buildtype=release',
  'fortran_std=legacy'
  ],
)

c_compiler = meson.get_compiler('c')

# BLAS
if host_machine.system() == 'darwin'
  message('Using `Accelerate` framework as BLAS library.')
  add_project_link_arguments('-framework', 'Accelerate', language: ['c', 'fortran'])
  blas_dep = declare_dependency()
else
  # static linking does not work in cibuildwheel
  blas_dep = dependency('openblas', required: false)
  if not blas_dep.found()
    blas_dep = c_compiler.find_library('openblas', required: false)
  endif
endif

if blas_dep.found()
  blas_source = []
else
  message('No specific BLAS library detected. Falling back to reference BLAS.')
  blas_source = join_paths(meson.current_source_dir(), 'extern/dblas.f')
endif

# Sources
daskr_source = files(
  'src/daskr_kinds.F90',
  'src/daskr_banpre.f90',
  'src/daskr_ilupre.f90',
  'src/daskr_rbdpre.f90',
  'src/daskr_rbgpre.f90',
  'src/daux.f',
  'src/daskr.f',
  'src/daskr_new.f90',
  'src/linpack.f90',
  'src/dsparskit.f',
  'src/blas_interfaces.f90',
)

daskr_lib = static_library(
  meson.project_name(),
  sources: [daskr_source, blas_source],
  dependencies: [blas_dep]
)

# incdir = include_directories('c/include')

# odrpack_dep = declare_dependency(
#   include_directories: [incdir],
#   link_with: [odrpack_lib]
# )

# Examples
if get_option('build_examples')
    subdir('example')
endif

# Tests
if get_option('build_tests')
  subdir('test')
endif