.. _contributing:

============
Contributing
============

This project is a community effort, and everyone is welcome to
contribute. Everyone within the community
is expected to abide by our `code of conduct <code_of_conduct.md>`_.

The project is hosted on
https://github.com/SteveDoyle2/pyNastran

Contributor Incubator
=====================

If you are interested in becoming a regular contributor to pyNastran, but
don't know where to start or feel insecure about it, you can join our non-public
communication channel for new contributors. To do so, please go to
https://discord.gg/s8RfSkZDHA and ask to be added to '#pyNastran'.
This is a general open source Discord channel for a variety of open source and commercial programs
and you can get guidance and support for your first few PRs.  This is a place you can ask questions
about anything: how to use git, github, how our PR review process works, technical questions
about the code, what makes for good documentation or a blog post, how to get involved involved
in community work, or get "pre-review" on your PR.


.. _submitting-a-bug-report:

Submitting a bug report
=======================

If you find a bug in the code or documentation, do not hesitate to submit a
ticket to the
`Issue Tracker <https://github.com/SteveDoyle2/pyNastran/issues>`_. You are
also welcome to post feature requests or pull requests.

If you are reporting a bug, please do your best to include the following:

1. A short, top-level summary of the bug. In most cases, this should be 1-2
   sentences.

2. A short, self-contained code snippet to reproduce the bug, ideally allowing
   a simple copy and paste to reproduce. Please do your best to reduce the code
   snippet to the minimum required.

3. The actual outcome of the code snippet.

4. The expected outcome of the code snippet.

5. The pyNastran version, Python version and platform that you are using. You
   can grab the version with the following commands::

      >>> import pyNastran
      >>> pyNastran.__version__
      '1.3.3'
      >>> import platform
      >>> platform.python_version()
      '3.7.7'

Thank you for your help in keeping bug reports complete, targeted and descriptive.

Requesting a new feature
========================

Please post feature requests to the
`Issue Tracker <https://github.com/SteveDoyle2/pyNastran/issues>`_.

The pyNastran developers will give feedback on the feature proposal. Since
pyNastran is an open source project with limited resources, we encourage
users to then also
:ref:`participate in the implementation <contributing-code>`.

.. _contributing-code:

Contributing code
=================

.. _how-to-contribute:

How to contribute
-----------------

The preferred way to contribute to pyNastran is to fork the `main
repository <https://github.com/SteveDoyle2/pyNastran/issues/>`__ on GitHub,
then submit a "pull request" (PR).

The best practices for using GitHub to make PRs to pyNastran are
documented in the :ref:`development-workflow` section.

A brief overview is:

1. `Create an account <https://github.com/join>`_ on GitHub if you do not
   already have one.

2. Fork the `project repository <https://github.com/SteveDoyle2/pyNastran>`_:
   click on the 'Fork' button near the top of the page. This creates a copy of
   the code under your account on the GitHub server.

3. Clone this copy to your local disk::

      $ git clone https://github.com/YourLogin/pyNastran.git

4. Create a branch to hold your changes::

      $ git checkout -b my-feature origin/main

   and start making changes. Never work in the ``main`` branch!

5. Work on this copy, on your computer, using Git to do the version control.
   When you're done editing e.g., ``pyNastran/bdf/bdf.py``, do::

      $ git add pyNastran/bdf/bdf.py
      $ git commit

   to record your changes in Git, then push them to GitHub with::

      $ git push -u origin my-feature

Finally, go to the web page of your fork of the pyNasran repo, and click
'Pull request' to send your changes to the maintainers for review.  You may
want to consider sending an email to the mailing list for more visibility.

.. seealso::

  * `Git documentation <https://git-scm.com/documentation>`_
  * `Git-Contributing to a Project <https://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project>`_
  * `Introduction to GitHub  <https://lab.github.com/githubtraining/introduction-to-github>`_
  * :ref:`development-workflow`
  * :ref:`using-git`

Contributing pull requests
--------------------------

It is recommended to check that your contribution complies with the following
rules before submitting a pull request:

* If your pull request addresses an issue, please use the title to describe the
  issue and mention the issue number in the pull request description to ensure
  that a link is created to the original issue.

* All public methods should have informative docstrings with sample usage when
  appropriate. Use the `numpy docstring standard
  <https://numpydoc.readthedocs.io/en/latest/format.html>`_.

In addition, you can should run the tests to check for errors with the following tools:

   * python all_tests.py

.. seealso::

  * :ref:`coding_guidelines`

.. _contributing_documentation:

Contributing documentation
==========================

You as an end-user of pyNastran can make a valuable contribution because you
more clearly see the potential for improvement than a core developer. For example, you can:

- Fix a typo
- Clarify a docstring
- Write or update an example
- Write or update a comprehensive tutorial

The documentation source files live in the same GitHub repository as the code.
Contributions are proposed and accepted through the pull request process.

For details see :ref:`how-to-contribute`.

If you have trouble getting started, you may instead open an `issue`_
describing the intended improvement.

.. _issue: https://github.com/SteveDoyle2/pyNastran/issues

.. _coding_guidelines:

Coding guidelines
=================

Adding new API
--------------

Every new function, parameter and attribute that is not explicitly marked as
private (i.e., starts with an underscore) becomes part of pyNastran's public
API. Changing the existing API is cumbersome. Therefore, take particular care when adding new API:

- Mark helper functions and internal attributes as private by prefixing them
  with an underscore.
- Carefully think about good names for your functions and variables.
- Try to adopt patterns and naming conventions from existing parts of the
  pyNastran API.

  __ https://emptysqua.re/blog/api-evolution-the-right-way/#adding-parameters


Thanks to https://github.com/matplotlib/matplotlib for their contributor's guide.  It's very well done and a great template!
