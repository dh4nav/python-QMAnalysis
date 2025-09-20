
YAML Format Reference
=====================

.. _yaml-format-reference:

This section describes the structure, options, and usage of the YAML input files for QMAnalysis. The format is inspired by the Docker Compose file reference and provides a comprehensive specification for each field.

Overview
--------

A QMAnalysis YAML file defines the workflow, input files, measurements, substitutions, and output. The following sections describe each top-level key, their options, and usage notes.

Top-level keys
--------------

.. list-table:: Top-level YAML keys
   :header-rows: 1

   * - Key
     - Type
     - Required
     - Description
   * - ``name``
     - string
     - yes
     - Name of the analysis procedure.
   * - ``comment``
     - string
     - no
     - Free-form comment or description.
   * - ``version``
     - integer
     - yes
     - Format version number. Used for compatibility.
   * - ``files``
     - list
     - yes
     - List of input file definitions. See :ref:`yaml-files-section`.
   * - ``measurements``
     - mapping
     - yes
     - Measurement definitions. See :ref:`yaml-measurements-section`.
   * - ``substitutions``
     - list
     - no
     - Substitution sets for atom replacements. See :ref:`yaml-substitutions-section`.
   * - ``output``
     - list
     - yes
     - Output actions (file, graph, etc). See :ref:`yaml-output-section`.

.. _yaml-files-section:

files
-----

Defines the input files for the analysis. Each entry is a mapping with the following options:

.. list-table:: ``files`` entry options
   :header-rows: 1

   * - Key
     - Type
     - Required
     - Description
   * - ``path``
     - string
     - yes
     - Path to the input file (relative or absolute).
   * - ``type``
     - string
     - yes
     - File type. Supported: ``xyz``, ``global_constants_csv``, ``gaussian_out``, etc.
   * - ``name``
     - string
     - no
     - Logical name for referencing this file elsewhere in the YAML.
   * - ``timestep``
     - string
     - no
     - Timestep label for trajectory files (e.g., ``init``, ``final``).

**Usage notes:**

- The ``name`` field is recommended if you have multiple files of the same type.
- The ``timestep`` field is used for time-dependent data (e.g., MD trajectories).

.. _yaml-measurements-section:

measurements
------------

Defines the measurements to compute. Each measurement type (e.g., ``distance``, ``angle``, ``dihedral``) is a key mapping to a list of measurement definitions.

.. list-table:: Supported measurement types
   :header-rows: 1

   * - Type
     - Description
   * - ``distance``
     - Pairwise distance between two atoms.
   * - ``angle``
     - Angle defined by three atoms.
   * - ``dihedral``
     - Dihedral angle defined by four atoms.

Each measurement entry supports:

.. list-table:: Measurement entry options
   :header-rows: 1

   * - Key
     - Type
     - Required
     - Description
   * - ``name``
     - string
     - yes
     - Name of the measurement (used for output and plotting).
   * - ``a``, ``b``, ``c``, ``d``
     - int or string
     - yes (as needed)
     - Atom indices (1-based) or aliases. Use as many as required for the measurement type.
   * - ``description``
     - string
     - no
     - Description of the measurement.

**Usage notes:**

- Atom indices are 1-based (first atom is 1).
- You may use aliases instead of indices if defined in the input files.
- The ``name`` field must be unique within each measurement type.

.. _yaml-substitutions-section:

substitutions
-------------

Defines sets of atom substitutions for systematic studies (e.g., mutating a residue or atom).

.. list-table:: ``substitutions`` entry options
   :header-rows: 1

   * - Key
     - Type
     - Required
     - Description
   * - ``name``
     - string
     - yes
     - Name of the substitution set.
   * - ``entries``
     - list
     - yes
     - List of substitution entries. Each entry must specify:
         - ``file`` (string, required): File name as defined in ``files``.
         - ``atom_index`` (int, required): 1-based atom index to substitute.

**Usage notes:**

- Substitutions are applied in the order listed.
- Useful for scanning over different atom replacements in a structure.

.. _yaml-output-section:

output
------

Defines the output actions to perform after analysis. Each entry is a mapping with one of the following keys:

- ``file``: Write results to a file (e.g., CSV).
- ``graph``: Generate a plot or graph.


.. list-table:: ``output`` entry options
   :header-rows: 1

   * - Key
     - Type
     - Required
     - Description
   * - ``file``
     - mapping
     - no
     - File output options:
         - ``path`` (string, required): Output file path.
         - ``type`` (string, required): Output type (e.g., ``csv``).
   * - ``graph``
     - mapping
     - no
     - Graph output options:
         - ``type`` (string, required): Graph type (e.g., ``scatter_plot``).
         - ``x``, ``y`` (string, required): Measurement names for axes.
         - ``file`` (string, required): Output file path.
         - ``file_format`` (string, optional): File format (e.g., ``tiff``, ``png``).
         - ``dpi`` (int, optional): Resolution in dots per inch.
         - ``font_family`` (string, optional): Font family for all text (e.g., ``Arial``, ``Times New Roman``).
         - ``font_size`` (int, optional): Base font size for labels and titles.
         - ``font_weight`` (string, optional): Font weight (e.g., ``normal``, ``bold``).

**Usage notes:**

- You may specify multiple output actions.
- Graph output supports various formats, DPI, and font settings for publication-quality figures.

Example
-------

.. code-block:: yaml

    name: Example Analysis
    comment: Demonstrates all YAML options
    version: 1
    files:
      - path: "constants.csv"
        type: global_constants_csv
      - path: "mol.xyz"
        type: xyz
        name: testmol
        timestep: "init"
    measurements:
      distance:
        - name: bond1
          a: 1
          b: 2
          description: "Test bond"
      angle:
        - name: angle1
          a: 1
          b: 2
          c: 3
    substitutions:
      - name: S1
        entries:
          - file: testmol
            atom_index: 1
    output:
      - file:
          path: ./out.csv
          type: csv
      - graph:
          type: scatter_plot
          x: bond1
          y: angle1
          file: ./plot.tiff
          file_format: tiff
          dpi: 300


See also
--------

- For API documentation, see the ``modules`` section in the navigation.
- For usage examples, see the project README or example YAML files in the repository.
