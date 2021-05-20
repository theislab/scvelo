Contributing guide
==================


Getting started
^^^^^^^^^^^^^^^

Contributing to scVelo requires a developer installation. As a first step, we suggest creating a new environment

.. code:: bash

    conda create -n ENV_NAME python=PYTHON_VERSION && conda activate ENV_NAME


Following, fork the scVelo repo on GitHub `here <https://github.com/theislab/scvelo>`.
If you are unsure on how to do so, please checkout the corresponding
`GitHub docs <https://docs.github.com/en/github/getting-started-with-github/quickstart/fork-a-repo>`.
You can now clone your fork of scVelo and install the development mode

.. code:: bash

    git clone https://github.com/YOUR-USER-NAME/scvelo.git
    cd scvelo
    git checkout --track origin/develop
    pip install -e '.[dev]'

The last line can, alternatively, be replaced by

.. code:: bash

    pip install -r requirements-dev.txt


Finally, to make sure your code follows our code style guideline, install pre-commit:

.. code:: bash

    pre-commit install


Coding style
^^^^^^^^^^^^

Our code follows `black` and `flake8` coding style. Code formatting (`black`, `isort`) is automated through pre-commit hooks. In addition, we require that

- functions are fully type-annotated.
- variables referred to in an error/warning message or docstrings are enclosed in \`\`.


Testing
^^^^^^^

To run the implemented unit tests locally, simply run

.. code:: bash

    python -m pytest


Documentation
^^^^^^^^^^^^^

The docstrings of scVelo largely follow the `numpy`-style. New docstrings should

- include neither type hints nor return types.
- reference an argument within the same docstrings using \`\`.


Submitting pull requests
^^^^^^^^^^^^^^^^^^^^^^^^

New features and bug fixes are added to the code base through a pull request (PR). To implement a feature or bug fix, create a branch from `develop`. For hotfixes use `master` as base. The existence of bugs suggests insufficient test coverage. As such, bug fixes should, ideally, include a unit test or extend an existing one. Please ensure that

- branch names have the prefix `feat/`, `bug/` or `hotfix/`.
- your code follows the project conventions.
- newly added functions are unit tested.
- all tests pass locally.
- if there is no issue solved by the PR, create one outlining what you try to add/solve and reference it in the PR description.
