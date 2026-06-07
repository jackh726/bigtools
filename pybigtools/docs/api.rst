Python API Reference
====================

.. currentmodule:: pybigtools

Opening files
-------------

.. autofunction:: open

Reading
-------

.. autoclass:: BBIReader
   :members:
   :member-order: groupwise

Writing
-------

.. autoclass:: BigWigWriter
   :members:

.. autoclass:: BigBedWriter
   :members:

Iterators
---------

.. autoclass:: BigWigIntervalIterator

.. autoclass:: BigBedEntriesIterator

Summary statistics
------------------

.. autoclass:: SummaryStatistics
   :members:

Exceptions
----------

.. autoexception:: BBIFileClosed

.. autoexception:: BBIReadError
