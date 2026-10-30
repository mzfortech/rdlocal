API reference
=============

All result objects are frozen dataclasses and provide ``to_dict()`` for JSON-friendly serialization. Rows with
non-finite input values are dropped. The treatment side is ``x >= cutoff``.

Estimation
----------

.. autofunction:: rdlocal.rd_estimate

``h`` and ``b`` accept a common value or a ``(left, right)`` pair. If ``h`` is omitted, ``bwselect`` chooses
the bandwidths. The default polynomial orders are ``p=1`` and ``q=p+1``. The ``fuzzy`` argument supplies
treatment uptake for a fuzzy design. ``kernel`` is ``"triangular"``, ``"epanechnikov"``, or ``"uniform"``;
``vce`` is ``"hc0"``–``"hc3"`` or ``"nn"``. Supplying ``cluster`` requests cluster-robust variance and makes
``vce`` irrelevant. ``level`` is the confidence level, such as ``0.95``.

.. autoclass:: rdlocal.RDResult
   :members:

Bandwidth selection
-------------------

.. autofunction:: rdlocal.rd_bandwidth

Supported selectors are ``mserd``, ``msetwo``, ``msesum``, ``msecomb1``, ``msecomb2``, ``cerrd``, ``certwo``,
``cersum``, ``cercomb1``, ``cercomb2``, and ``ik``. The IK selector is implemented for sharp designs with
``p=1`` and the triangular kernel. ``p`` is the estimation order, ``q`` is the higher bias-estimation order,
``kernel`` selects the kernel, ``cluster`` supplies optional cluster identifiers, and ``scaleregul`` sets the
bandwidth selector's regularization factor.

.. autoclass:: rdlocal.BandwidthResult
   :members:

Density manipulation test
--------------------------

.. autofunction:: rdlocal.rd_density_test

``h`` may be omitted for a data-driven common bandwidth, or set to a common value or ``(left, right)`` pair.
``p`` is the density polynomial order; ``kernel`` is ``"triangular"``, ``"epanechnikov"``, or ``"uniform"``.

.. autoclass:: rdlocal.DensityResult
   :members:

Binned means
------------

.. autofunction:: rdlocal.rd_binned

``nbins`` can be one integer for both sides or a ``(left, right)`` pair. ``spacing`` is ``"evenly"`` or
``"quantile"``. The cutoff separates the two sides; non-finite input rows are omitted.

.. autoclass:: rdlocal.BinnedResult
   :members:
