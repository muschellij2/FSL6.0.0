0.2.1 (Friday 12th October 2018)
--------------------------------


* Adjustments to ``AWBase::operator<<`` operator overloads.


0.2.0 (Thursday 4th October 2018)
---------------------------------


* Fixes to ``.Storage`` calculation on certain subviews
* Refactored tests, and added a ``run_tests.sh`` script to make testing
  easier. Added GitLab CI integration for automated testing on pushes.
* Newmat source code is now included in the ``tests`` directory.


0.1.0 (Tuesday 2nd October 2018)
--------------------------------


* Fixed bug in subview lookup/assignment for non-simple matrix types
* Added ability to call ``.Storage`` on sub-views (except for band matrix
  types).


0.0.0
-----

* Initial working implementation
