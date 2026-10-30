Getting started
===============

Install
-------

``rdlocal`` requires Python 3.11 or newer:

.. code-block:: console

   python -m pip install rdlocal

Quick example
-------------

Estimate the jump in an outcome at a cutoff. By default, observations with a running-variable value at or above
the cutoff are treated, and the function selects bandwidths automatically.

.. code-block:: python

   import numpy as np
   from rdlocal import rd_estimate

   rng = np.random.default_rng(0)
   x = rng.uniform(-1, 1, 2000)
   y = 0.5 * x + 0.4 * (x >= 0) + rng.normal(0, 0.3, 2000)

   result = rd_estimate(y, x, cutoff=0.0)
   print(result.coef, result.coef_bc, result.ci_robust)

``coef`` is the conventional estimate; ``coef_bc`` and ``ci_robust`` report the bias-corrected estimate and
robust confidence interval. See the :doc:`API reference <api>` for manual bandwidths, fuzzy designs, and other
options.
