fenicsR13
=========

.. image:: media/logo.svg
    :width: 200px
    :alt: logo
    :align: left

.. inclusion-marker

**fenicsR13**: *A Tensorial Mixed Finite Element Solver for the Linear R13 Equations Using the FEniCS Computing Platform*

|pipeline| |coverage| |version| |docs_website| |zenodo|

.. list-table::
   :widths: 20 50
   :header-rows: 0

   * - Paper
     - **Lambert Theisen**, **Manuel Torrilhon**. 2021. *fenicsR13: A Tensorial Mixed Finite Element Solver for theLinear R13 Equations Using the FEniCS Computing Platform*. ACM Trans. Math. Softw. 47, 2, Article 17 (April 2021), 29 pages, `10.1145/3442378 <https://dx.doi.org/10.1145/3442378>`_. `Preprint @ arXiv (Free PDF) <https://dl.acm.org/doi/10.1145/3442378?cid=99659694875>`_.

Main Features
--------------------------------------------------------------------------------

- Solving Steady 2D Linear R13 Equations in Normalized Form
- Stabilization with: Continous Interior Penalty (CIP) or Galerkin Least Squares (GLS)
- Option for convergence study with input for exact solution
- Output in XDMF_/HDF5 and PDF
- Arbitray finite element combinations (thanks to FEniCS_)
- Easy setup of parameter studies
- Structured/Documented YAML_ input file format
- Rarefied gas flow effects predictable: Knudsen paradox, Knudsen pump, ...
- Interface for Gmsh_ meshes
- Full exposure of the PETSc_ preconditioned iterative solver library

.. _FEniCS: https://fenicsproject.org/
.. _YAML: https://de.wikipedia.org/wiki/YAML
.. _XDMF: http://www.xdmf.org/index.php/XDMF_Model_and_Format
.. _Gmsh: http://gmsh.info/
.. _PETSc: https://petsc.org/

Installation
--------------------------------------------------------------------------------

1. Install `Git <https://git-scm.com/>`_ (if not shipped with you OS).
2. Install `Git LFS <https://git-lfs.com/>`_ (needed to download the meshes from the repository).
3. Clone and open the repository:

.. code-block:: bash

    # Clone Repository and open main folder
    git clone https://git.rwth-aachen.de/lamBOO/fenicsR13.git
    cd fenicsR13

*Alternatively*: `download <https://git.rwth-aachen.de/lamBOO/fenicsR13/-/archive/master/fenicsR13-master.zip>`_ the repository, unpack and open it.

Docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is recommended to use the program with `Docker`_, a virtualization software available (with free docker account) for `download`_ for most operating systems. Our repository provides a Docker container based on the official FEniCS Docker image which contains all necessities for FEniCS and has been adjusted to solve the R13 equations. The program `docker-compose`_ is also required.

.. _`Docker`: https://en.wikipedia.org/wiki/Docker_(software)
.. _`download`: https://www.docker.com/products/docker-desktop
.. _`docker-compose`: https://docs.docker.com/compose/install/

Starting the Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure you installed Docker and it is running on your system. You can start the fenicsR13 compute environment from the root-directory of the repository using the command-line

.. code-block:: bash

    ## 1) FOR USERS:
    ## - Pulls and runs "fenicsr13_release" service
    ## - Uses preinstalled latest version of fenicsR13
    docker compose pull fenicsr13_release
    docker compose run --rm fenicsr13_release

    ## 2) FOR DEVELOPERS (uncommend the last 3 lines of this block):
    ## - Builds the "fenicsr13_debug"
    ## - Install an editable version of fenicsR13
    ## - The "pip"-command install an editable version of fenicsR13
    ## - This allows a modification of the source files which can be directly
    ##   executed within the debug container
    # docker compose build fenicsr13_debug
    # docker compose run fenicsr13_debug
    # sudo pip install -e .

When you run this for the first time, docker will pull (download and extract) the container image from our repository which is roughly 800MB and the download may require some patience. After the initial download the docker image will be stored (2-3 GB) on your system and any new run will start the container immediately.

The above docker command will start a RedHat-based bash environment. The repository is mounted as a volume under ``/home/fenics/shared`` in the container and should be the default folder on startup. What ever is changed in this folder or its sub-folders will change in the folder on your original system as well. However, any changes outside this folder will be gone after you exit the docker environment.

Running a Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To execute a simulation case, go to to the case folder (e.g. ``examples/lid_driven_cavity``)

.. code-block:: bash

    # [It might be required to install the package inside the Docker container]
    pip install --user -e .

.. code-block:: bash

    # Move to folder:
    cd examples/lid_driven_cavity

We provide a script to utilize `gmsh`_ and generate a `H5`_ mesh-file from a local geometry file by

.. code-block:: bash

    # Create mesh:
    ./create_mesh.sh

To run a simulation execute the solver main program ``fenicsR13.py`` (which is located in the ``src``-directory in the top level) while specifying an input file as first command line argument.

.. code-block:: bash

    # Run program with given input file:
    fenicsR13 input.yml


Output files will be written to a folder which is named after the ``output_folder`` keyword of the ``input.yml``. For immediate inspection the output folder contains simple visualizations in PDF files for each of the fields (temperature, pressure,...).

The numerical results for each field is ouput into ``h5``-files, including mesh data and with corresponding ``xdmf``-file. The XDMF-files can be opened in Paraview to perform visualization, e.g., with ``Paraview > File > Open > u_0.xdmf > Apply filters``

.. _`gmsh`: http://gmsh.info/
.. _`H5`: https://en.wikipedia.org/wiki/Hierarchical_Data_Format

.. code-block:: bash

    # Leave directory:
    cd ../..

**Channel Flow Example**

We provide a simple example of a flow through a finite-length channel in 2D.

.. code-block:: bash

    # Move to folder:
    cd examples/channel_flow_force
    # Create mesh:
    ./create_mesh.sh
    # Run program with given input file:
    fenicsR13 input.yml

In the output folder the results can be post-processed to demonstrate the `Knudsen paradox`_ in a simple table.

.. code-block:: bash

    # Go to folder with simulation results (=casename in input.yml)
    cd channel_flow_force
    # Generate correlation data between Knudsen number and massflow
    bash postprocessing.sh
    cat table.csv
    # Leave directory:
    cd ../..

.. _`Knudsen paradox`: https://en.wikipedia.org/wiki/Knudsen_paradox

**Convergence Study**

We can test the convergence of the R13 discretization on a simple double-cylindrical geometry.

.. code-block:: bash

    # Move to folder:
    cd tests/2d_r13
    # Meshes are already in Git:
    ls ../2d_mesh
    # Run program with given input file:
    fenicsR13 inputs/r13_1_coeffs_nosources_norot_inflow_p1p1p1p1p1_stab.yml
    # Go to folder with simulation results (=casename in input.yml)
    cd r13_1_coeffs_nosources_norot_inflow_p1p1p1p1p1_stab
    # Open errors:
    cat errors.csv



Additional information
--------------------------------------------------------------------------------

Parallel Execution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FEniCS allows simple parallelization using MPI

.. code-block:: bash

    # Parallel execution ("-u" to flash stdout)
    # Usage: mpirun -n <numberOfProcesses> <serialCommand>
    # E.g.: mpirun -n 4 fenicsR13 input.yml

Building the Docker Image Locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main folder of this repository contains a ``Dockerfile`` defining the used environment. Here, we used the optimized and official FEniCS Docker image and include ``Gmsh`` and install some requirements from the ``requirements.txt``. This can take a while, especially the ``Gmsh`` mirror can be quite slow. To avoid very long execution commands (``docker run <..> -v <volume share> <etc..>``), a ``docker-compose.yml`` is used to store all these parameters. ``docker compose`` acts as an wrapper for the Docker execution.

The ``fenics`` environment (also called *service* in the ``docker-compose.yml``) first has to be build and can be executed afterwards. The command to build the container is

.. code-block:: bash

    # build fenics service
    docker compose build fenicsr13_release


Interactive Docker Sessions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to use a Jupyter sever or a X11 forwarding but this is not recommended anymore. All relevant plots are now written by default without the need for the tricky X11 forwarding or interactive usage with Jupyter.

Documentation
--------------------------------------------------------------------------------

Documentation using Sphinx is available.

Pre-Build Version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Visit the hosted version on `Gitlab Pages`_ or download the artifacts from Gitlab's CI ``pages``-pipeline.

.. _`Gitlab Pages`: https://lamboo.pages.rwth-aachen.de/fenicsR13/

Manual Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # cat .gitlab-ci.yml
    cd docs
    sphinx-apidoc -o source/src ../src
    sphinx-apidoc -o source/tests/2d_heat ../tests/2d_heat
    sphinx-apidoc -o source/tests/2d_stress ../tests/2d_stress
    sphinx-apidoc -o source/tests/2d_r13 ../tests/2d_r13
    sphinx-apidoc -o source/tests/3d_heat ../tests/3d_heat
    sphinx-apidoc -o source/tests/3d_stress ../tests/3d_stress
    sphinx-apidoc -o source/tests/3d_r13 ../tests/3d_r13
    sphinx-apidoc -o source/examples ../examples
    make html
    make latex

Developer Legacy Notes
--------------------------------------------------------------------------------

Developer Tips
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Monitor the performance of the program with e.g.:

    .. code-block:: bash

        htop -p `{ fenicsR13 inputs/1_coeffs_nosources_norot_inflow_p1p1p1p1p1_stab.yml > /dev/null & } && echo $!`

- Use doctest with ``python3 -m doctest -v src/meshes.py``
- Run ``pydocstyle`` once in a while
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

Python notes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Get current work directory:

    .. code-block:: python

        import os
        cwd = os.getcwd()
        print(cwd)

- Latex font for matplotlib:

    .. code-block:: python

        # LaTeX text fonts:
        # Use with raw strings: r"$\mathcal{O}(h^1)$"
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

- Get system path where modules are searched:

    .. code-block:: python

        import sys
        print(sys.path)

Create new version tag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Add CHANGELOG entry
2. Adapt version in `conf.py` for docs and `setup.py` for package
3. Change badge in ``README.rst``
4. Change version in program information printing
5. Build new Docker container

Gitlab CI Setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- The ``build`` stage has to be triggered manually when something in the setup changes. This is because it takes a fair amount of time.
- In ``~/.gitlab-runner/config.toml`` (for the runner):
    - change priviliges to true
    - Use local images: ``pull_policy = "if-not-present"``
    - To ``[[runners]]`` add ``environment = ["DOCKER_TLS_CERTDIR="]`` (See https://gitlab.com/gitlab-org/gitlab-ce/issues/64959)
- Run local: ``gitlab-runner exec docker --docker-privileged build`` or with ``build`` replaced by job name
    - maybe local vars have to be change to use local Docker images because ``CI_REGISTRY``,... are not set

An example gitlab runner ``config/toml`` in ``~/.gitlab-runner`` can look like:

.. code-block:: toml

    concurrent = 1
    check_interval = 0

    [[runners]]
    name = "190716-macbookpro"
    url = "https://git.rwth-aachen.de/"
    token = "<PRIVATE_TOKEN>"
    executor = "docker"
    environment = ["DOCKER_TLS_CERTDIR="]
    [runners.docker]
        tls_verify = false
        image = "docker:stable"
        privileged = true
        disable_cache = false
        volumes = ["/cache"]
        shm_size = 0
        pull_policy = "if-not-present"
    [runners.cache]

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

This is may be a convenient solution.
Run a file with ``%run ../../fenicsr13/fenicsr13.py``

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

Contact
--------------------------------------------------------------------------------

:Author:
    | Lambert Theisen
    | lambert.theisen@rwth-aachen.de
:Supervisor:
    | Prof. Dr. Manuel Torrilhon
    | Lehrstuhl für Mathematik (MathCCES)
    | RWTH Aachen University
    | mt@mathcces.rwth-aachen.de

.. |pipeline| image:: https://git.rwth-aachen.de/lamboo/fenicsR13/badges/master/pipeline.svg
    :target: https://git.rwth-aachen.de/lamboo/fenicsR13/commits/master
    :alt: Pipeline status

.. |coverage| image:: https://git.rwth-aachen.de/lamboo/fenicsR13/badges/master/coverage.svg
    :target: https://git.rwth-aachen.de/lamboo/fenicsR13/pipelines
    :alt: Test coverage

.. |version| image:: https://img.shields.io/badge/version-1.5-blue.svg
    :target: https://git.rwth-aachen.de/lamBOO/fenicsR13/-/tags
    :alt: Documentation Website

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3673039.svg
    :target: https://doi.org/10.5281/zenodo.3673039
    :alt: Zenodo Link with DOI

.. |docs_website| image:: https://img.shields.io/badge/docs-dev-blue.svg
    :target: https://lamboo.pages.rwth-aachen.de/fenicsR13/
    :alt: Documentation Website
