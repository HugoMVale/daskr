examples = [
    'heat',
    'heatilu',
    'web'
    # 'webilu'
]

foreach example : examples
    executable(
      'example_@0@'.format(example),
      sources: files('example_@0@.f90'.format(example.underscorify())),
      link_with: daskr_lib,
    )
endforeach