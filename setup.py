from setuptools import setup

setup(
    name='vcf_from_snvphyl',
    version='0.0.1',
    author='Peter van Heusden',
    author_email='pvh@sanbi.ac.za',
    url='https://github.com/COMBAT-TB/vcf_from_snvphyl',
    project_urls={"Bug Tracker": 'https://github.com/COMBAT-TB/vcf_from_snvphyl/issues'},
    description='Parses variant report from SNVPhyl and generates annotated (ala. SnpEff) VCF files',
    keywords='neo4j,COMBAT-TB,vcf,tb',
    license="MIT",
    packages=['snptools'],
    setup_requires=[
        'pytest-runner',
    ],
    tests_require=[
        'pytest',
    ],
    install_requires=[
        'biopython',
        'intervaltree',
        'py2neo',
        'pandas',
        'tqdm',
    ],
    scripts=[
        'bin/vcf_from_snvphyl',
    ],
)
