from setuptools import setup

setup(
    name='EukMetaSanity',
    version='0.1',
    description='Eukaryotic Metagenome-Assembled Genome Structural/Functional Annotation Pipeline',
    url='https://github.com/cjneely10/EukMetaSanity',
    author='Christopher Neely',
    author_email='christopher.neely1200@gmail.com',
    license='GPL-3.0',
    packages=['EukMetaSanity'],
    scripts=["EukMetaSanity/EukMetaSanity"],
    zip_safe=False
)
