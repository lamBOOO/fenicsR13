.. image:: logo.png
    :width: 200px
    :alt: logo

.. inclusion-marker

gasDynamicsFEM
================================================================================

|pipeline| |coverage|

``#extendedGasDynamics`` ``#using`` ``#FEniCS``

  Repository for Master thesis project regarding FEM simulations for
  non-equilibrium gas dynamics.

Installation
--------------------------------------------------------------------------------

It is recommended to use the program within a Docker container.

Docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install `Docker Desktop`_ for your OS.

.. _`Docker Desktop`: https://www.docker.com/products/docker-desktop

.. code-block:: console

    # build and run fenics service
    docker-compose build fenics
    docker-compose run fenics

    # Execute some examples
    cd fenicsR13/examples/heat
    python3 ../../src/fenicsR13.py


The main folder of this repository contains a ``Dockerfile`` defining the used environment. Here, we used the optimized and official FEniCS Docker image and include ``Gmsh`` and install some requirements from the ``requirements.txt``. This can take a while, especially the ``Gmsh`` mirror can be quite slow. To avoid very long execution commands (``docker run <..> -v <volume share> <etc..>``), a ``docker-compose.yml`` is used to store all these parameters. ``docker-compose`` acts as an wrapper for the Docker execution.

The ``fenics`` environment (also called *service* in the ``docker-compose.yml``) first has to be build and can be executed afterwards. The steps to perform then read

The whole repository is mounted as a volume under ``/home/fenics/shared`` in the container and should be the default folder on startup. To execute the solver, move to the case folder (e.g. ``/home/fenics/shared/cases/heatSystem``) and execute the script (e.g. ``python3 heat.py``). Output files should be written in that case, e.g. to the ``results`` folder.

It is convenient to use a Jupyter sever or a X11 forwarding.

macOS Native FEniCS Installation (not recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Install ``miniconda`` from `here <https://conda.io/projects/conda/en/latest/user-guide/install/macos.html>`_
    #. If using ``zsh``, add miniconda bins to PATH: ``export PATH="$HOME/ miniconda3/bin:$PATH"`` to ``~/.zshrc``
    #. Maybe, activation has to be done with executing ``<path to  miniconda>/bin/activate``
    #. Optional: Create separate coda environment: ``conda creafenics-env``
#. Install FEniCS using conda: ``conda install -c conda-forge fenics``
    #. Optional: Install ``matplobib``: ``conda install -c conda-forge  matplotlib``
    #. Optional: Install ``meshio``: ``conda install -c mrossi meshio``
    #. Optional (for linting): ``conda install pylint``
    #. Install mshr with ``conda install -c conda-forge mshr``
    #. Fix macOS bug in matplotbib: ``mkdir -p ~/.matplotlib; echo  "backend: TkAgg" > ~/.matplotlib/matplotlibrc``
    #. XCode and command line developer tools msut be installed!
    #. Optional: Install Jupyter: ``conda install -c anaconda jupyter``
    #. Optional: Install documentation system: ``conda install -c anaconda  sphinx``
    #. Optional: ``conda install -c anaconda sympy``

Further Installation Tips
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Interactive Jupyter Notebooks with Microsoft's Visual Studio Code**

This is the most convenient solution.
Run a file with ``%run ../../src/fenicsr13.py``

**X11 Window Forwarding on OSX**

See guide_ for the programs to install. Then source the ``open-macos-gui-tunnel.sh`` with ``. open-macos-gui-tunnel``. Afterwards, start the container and run the ``change-matplotbib-backend-tkagg.sh`` script to set the right ``matplotlib``'s output.

.. _guide: http://joshuamccall.com/articles/docker.html

**X11 Window Forwarding on Windows**

This has to be studied.

Documentation
--------------------------------------------------------------------------------

Documentation using Sphinx is available.

Manual Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: console

    cd docs
    sphinx-apidoc -o source/ ../src
    make html
    # open _build/html/index.html

Pre-Build Version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download pre-created artifacts from Gitlab's CI pipeline page. Or visit the
hosted version on `Gitlab Pages`_.

.. note:: Currently, there's a bug regarding Gitlab Pages with internal repositories. The access control fails and the page cannot be accessed. This may be fixed in the future.

.. _`Gitlab Pages`: https://lamboo.pages.rwth-aachen.de/gasdynamicsfem/

Developer Notes
--------------------------------------------------------------------------------

- Matplotbib fails when having wrong backend on macOS
    - Fix: Add ``backend: TkAgg`` to ``~/.matplotlib/matplotlibrc`` file
- Performance in Docker is way better than conda build, especially JIT compilation is faster
- Get C++ inlcude paths: ``echo | gcc -E -Wp,-v -``
- Bessel functions in DOLFIN:
    - C++17 functions cannpot be used. Boost functions also not per default. ``Expression("boost::math::cyl_bessel_i(0,atan2(x[1], x[0]))", degree=2)`` is allowed if one changes in file ``/usr/local/lib/python3.6/dist-packages/dolfin/jit/jit.py``

        .. code-block:: python

            _math_header = """
            // cmath functions
            #include <boost/math/special_functions/bessel.hpp> // Added
            %s
            """

- Python notes:
    - Get current work directory

        .. code-block:: python

            import os
            cwd = os.getcwd()
            print(cwd)

    - Latex font for matplotlib

        .. code-block:: python

            # LaTeX text fonts:
            # Use with raw strings: r"$\mathcal{O}(h^1)$"
            # plt.rc('text', usetex=True)
            # plt.rc('font', family='serif')

    - Get system path where modules are searched

        .. code-block:: python

            import sys
            print(sys.path)

- Gitlab CI Setup:
    - In ``~/.gitlab-runner/config.toml`` (for the runner):
        - change priviliges to true
        - Use local images: ``pull_policy = "if-not-present"``
    - Run local: ``gitlab-runner exec docker --docker-privileged build`` or with ``build`` replaced by job name
        - maybe local vars have to be change to use local Docker images because ``CI_REGISTRY``,... are not set

Contact
-------

.. image:: media/mathcces.png
    :width: 400px
    :alt: logo

:Author:
    | Lambert Theisen
    | lambert.theisen@rwth-aachen.de
:Supervisor:
    | Prof. Dr. Manuel Torrilhon
    | Lehrstuhl fur Mathematik (MathCCES)
    | RWTH Aachen University
    | mt@mathcces.rwth-aachen.de
:Context:
    | Masterthesis Computational Engineering Science
    | RWTH Aachen University
    | *Simulation of Non-Equilibrium Gas Flows Using the FEniCS Computing Platform*

.. |pipeline| image:: https://git.rwth-aachen.de/lamboo/gasdynamicsfem/badges/master/pipeline.svg
    :target: https://git.rwth-aachen.de/lamboo/gasdynamicsfem/commits/master
    :alt: Pipeline status

.. |coverage| image:: https://git.rwth-aachen.de/lamboo/gasdynamicsfem/badges/master/coverage.svg
    :target: https://git.rwth-aachen.de/lamboo/gasdynamicsfem/commits/master
    :alt: Test coverage
