Limitations and validation
==========================

The estimators agree closely with published results and black-box reference values when evaluated at the same
bandwidths. Some data-driven bandwidth selectors still differ from reference implementations because parts of
their pilot procedures are not fully specified in the available papers. The current choices and measured gaps
are documented in the project's `clean-room record <https://github.com/mzfortech/rdlocal/blob/main/CLEANROOM.md>`_.

For results that need to be compared across software, supply ``h`` and ``b`` explicitly. The Senate vote-margin
data used for one published example is not distributed with the package.
