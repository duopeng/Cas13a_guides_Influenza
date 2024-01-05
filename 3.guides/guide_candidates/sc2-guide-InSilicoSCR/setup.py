from setuptools import setup

setup(name='sc2-guide-InSilicoSCR',
      version="0.0.1",
      summary='SC2 Guides InSilico Screening Tools',
      description='https://github.com/czbiohub/sc2-guide-InSilicoSCR/wiki/InSilicoSCR',
      url='https://github.com/czbiohub/sc2-guide-InSilicoSCR',
      author='Chunyu Zhao',
      author_email='chunyu.zhao@czbiohub.org',
      license='MIT',
      packages=['isscrlib'],
      install_requires=[],
      dependency_links=[],
      entry_points={
        'console_scripts': [
          'isscr_gen_neighbors = isscrlib.gen_neighbors:main'
        ]
      },
      zip_safe=False
)
