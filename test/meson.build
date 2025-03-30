tests = [
  'krdem1',
  'krdem2'
]

foreach t : tests
  test(
    t,
    executable(
      'test-@0@'.format(t),
      sources: files('test_@0@.f90'.format(t.underscorify())),
      link_with: daskr_lib,
    ),
    suite: 'unit',
  )
endforeach