.. overview:

API Overview
============

**EukMetaSanity** has a simple set of high-level goals:

* Run a series of analysis programs on a set of input files.
* Parallelize to match system resources.
* Collect all output data and package for user.

Each of these goals is accomplished through the ``TaskManager`` class:

.. autoclass:: EukMetaSanity.tasks.base.task_manager.TaskManager
	:members:

