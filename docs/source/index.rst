EukMetaSanity
=============
Highly parallelized platform for eukaryote annotation
-----------------------------------------------------

Biological Data Analysis Pipelines
==================================
**EukMetaSanity** is packaged with three pipelines:

The ``run`` pipeline is the first step in eukaryotic annotation. Specifically, it will:

* Provide an initial taxonomy assignment by searching combined OrthoDB and MMETSP databases.
* Mask repetitive regions using data from NCBI and/or ab-initio predicted repeats.
* Generate ab-initio predictions of exon/intron structure and boundaries.
* Integrate protein evidence from OrthoDB and MMETSP to validate or extend ab-initio predictions.

After initial evidence-based gene models are created, users may elect to identify functional annotations using
``report``:

* Search EggNOG database.
* Search KEGG database.
* Search OrthoDB/MMETSP databases.
* Search any mmseqs sequence or profile database.

Users may also elect to integrate transcriptomic evidence into their prediction pipeline using the ``refine`` step:

* Map RNA-seq data or assembled transcriptomes to genome/MAG.

Data Pipeline API
=================
**EukMetaSanity** is built upon a flexible API combining resource acquisition/task scheduling operations from Python's
``dask`` with filesystem access/program calling functions of ``plumbum``. This API is available for users who wish to
generate their own data analysis pipelines using **EukMetaSanity's** backend.

.. toctree::
	:maxdepth: 2

	api/overview
	api/pipeline_manager
	api/task
	api/example

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
