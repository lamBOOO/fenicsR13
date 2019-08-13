.. image:: logo.png
    :width: 200px
    :alt: logo

.. inclusion-marker

fenicsR13
================================================================================

|pipeline| |coverage| |version| |website|

``#extendedGasDynamics`` ``#using`` ``#FEniCS``

  Repository for Master thesis project regarding FEM simulations for
  non-equilibrium gas dynamics.

Installation
--------------------------------------------------------------------------------

It is recommended to use the program within a Docker container.

Docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First install `Docker Desktop`_ for your OS. Then:

.. _`Docker Desktop`: https://www.docker.com/products/docker-desktop

.. code-block:: bash

    # build and run fenics service
    docker-compose build fenics
    docker-compose run fenics

    # Execute some examples
    cd fenicsR13/examples/heat

    # in serial...
    # usage: python3 <path/to/fenicsR13.py> <path/to/input.yml>
    python3 ../../src/fenicsR13.py inputs/01_coeffs_p1p1_stab.yml

    # in parallel... ("-u" to flash stdout)
    # usage: mpirun -n <numberOfProcesses> <serialCommand>
    mpirun -n 4 python3 -u ../../src/fenicsR13.py inputs/01_coeffs_p1p1_stab.yml

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

A nice guide can be found `here on Dev.to`_.

.. _`here on Dev.to`: https://dev.to/darksmile92/run-gui-app-in-linux-docker-container-on-windows-host-4kde

The steps can be summarized as:

1. Install the package manager `Chocolatey`_.

    .. code-block:: dosbatch

        REM comment: open cmd.exe as admin
        @"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -InputFormat None -ExecutionPolicy Bypass -Command "iex ((New-Object System.Net.WebClient).DownloadString('https://chocolatey.org/install.ps1'))" && SET "PATH=%PATH%;%ALLUSERSPROFILE%\chocolatey\bin"

2. Open ``cmd.exe`` as admin and install `VcXsrv Windows X Server`_.

    .. code-block:: bash

        choco install vcxsrv
3. Open a X11 server and set the ``ip`` variable (that is used in the ``docker-compose.yml`` when starting the Docker container to set ``export DISPLAY=${ip}:0``).

    .. code-block:: bash

        # home of this repo
        source sripts/open-windows-gui-tunnel.sh

.. _`Chocolatey`: https://chocolatey.org/
.. _`VcXsrv Windows X Server`: https://sourceforge.net/projects/vcxsrv/

Documentation
--------------------------------------------------------------------------------

Documentation using Sphinx is available.

Manual Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    cd docs
    sphinx-apidoc -o source/ ../src
    make html
    # open _build/html/index.html

Pre-Build Version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download pre-created artifacts from Gitlab's CI pipeline page. Or visit the
hosted version on `Gitlab Pages`_.

.. note:: Currently, there's a bug regarding Gitlab Pages with internal repositories. The access control fails and the page cannot be accessed. This may be fixed in the future.

.. _`Gitlab Pages`: https://lamboo.pages.rwth-aachen.de/fenicsR13/

Developer Notes
--------------------------------------------------------------------------------

- Monitor the performance of the program with e.g.:

    .. code-block:: bash

        htop -p `{ python3 ../../src/fenicsR13.py inputs/1_coeffs_nosources_norot_inflow_p1p1p1p1p1_stab.yml > /dev/null & } && echo $!`

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
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')

    - Get system path where modules are searched

        .. code-block:: python

            import sys
            print(sys.path)

- Create new version tag:
    1. Add CHANGELOG entry
    2. Adapt version in `conf.py` for docs

- Gitlab CI Setup:
    - The ``build`` stage has to be triggered manually when something in the setup changes. This is because it takes a fair amount of time.
    - In ``~/.gitlab-runner/config.toml`` (for the runner):
        - change priviliges to true
        - Use local images: ``pull_policy = "if-not-present"``
        - To ``[[runners]]`` add ``environment = ["DOCKER_TLS_CERTDIR="]`` (See https://gitlab.com/gitlab-org/gitlab-ce/issues/64959)

    - Run local: ``gitlab-runner exec docker --docker-privileged build`` or with ``build`` replaced by job name
        - maybe local vars have to be change to use local Docker images because ``CI_REGISTRY``,... are not set

Contact
-------

.. image:: ./media/mathcces.png
    :width: 400px
    :alt: mathcces
    :target: http://www.mathcces.rwth-aachen.de

:Author:
    | Lambert Theisen
    | lambert.theisen@rwth-aachen.de
:Supervisor:
    | Prof. Dr. Manuel Torrilhon
    | Lehrstuhl für Mathematik (MathCCES)
    | RWTH Aachen University
    | mt@mathcces.rwth-aachen.de
:Context:
    | Masterthesis Computational Engineering Science
    | RWTH Aachen University
    | *Simulation of Non-Equilibrium Gas Flows Using the FEniCS Computing Platform*

.. |pipeline| image:: https://git.rwth-aachen.de/lamboo/fenicsR13/badges/master/pipeline.svg
    :target: https://git.rwth-aachen.de/lamboo/fenicsR13/commits/master
    :alt: Pipeline status

.. |coverage| image:: https://git.rwth-aachen.de/lamboo/fenicsR13/badges/master/coverage.svg
    :target: https://git.rwth-aachen.de/lamboo/fenicsR13/pipelines
    :alt: Test coverage

.. |version| image:: https://img.shields.io/badge/version-v0.3-blue.svg
    :target: https://git.rwth-aachen.de/lamBOO/fenicsR13/-/tags
    :alt: Documentation Website

.. |website| image:: https://img.shields.io/badge/doc-https%3A%2F%2Flamboo.pages.rwth--aachen.de%2FfenicsR13%2F-blue.svg
    :target: https://lamboo.pages.rwth-aachen.de/fenicsR13/
    :alt: Documentation Website
