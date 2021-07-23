from setuptools import setup

setup(
    name='EukMetaSanity',
    version='0.1.0',
    description='Eukaryotic Metagenome-Assembled Genome Structural/Functional Annotation Pipeline',
    url='https://github.com/cjneely10/EukMetaSanity',
    author='Christopher Neely',
    author_email='christopher.neely1200@gmail.com',
    license='GPL-3.0',
    packages=['EukMetaSanity', "EukMetaSanity.data"],
    scripts=["EukMetaSanity/download-data"],
    zip_safe=False
)
