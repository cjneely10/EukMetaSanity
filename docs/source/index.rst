EukMetaSanity: Highly parallelized platform for eukaryote annotation
====================================================================

Biological Data Analysis Pipelines
==================================
**EukMetaSanity** is packaged with three pipelines:

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   run
   report
   refine

``run`` accepts a directory of FASTA files and outputs a series of gene annotations in GFF3 and FASTA format.

Eukaryotic gene annotation requires additional calculation that is not present in prokaryotic annotation. Specificall,
the ``run`` pipeline will:

* Provide an initial taxonomy assignment

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
