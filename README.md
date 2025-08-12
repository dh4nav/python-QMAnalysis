# Input File Format

## Top-Level Elements

* name
* comment
* version
* files
  * 
    * path
    * name#
    * type
    * options#
* substitutions#
  * name
  * entries
    * file
    * atom_index
* sequences#
  * name
  * entries
    * identifier
    * ...#
* measurements#
  * distance
    * name
    * a
    * b
  * angle
    * name
    * a
    * b
    * c
  * dihedral
    * name
    * a
    * b
    * c
    * d
* graphics
  * plot
    * name
    * identifier
  * ...#