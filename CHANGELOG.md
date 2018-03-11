# v0.3 (10 Mar 2018)
* Contacts may now also be polygons.
* Diff, poly, and metal may now also be rectangles.
* The cubic Bezier curves c, C are now accepted in paths, and handled by approximating with three segments (four points).
* Self-crossing paths are detected and verbosely reported for correction.
* SVG multi-shell and multi-hole paths are now properly handled.
* n-input NOR gates and power NOR gates are now found (a 1-input NOR gate is an inverter)
    * The test/qnames.svg file contains an inverter.
    * Note that drivers formed from parallel transistors are not yet handled.
* Any signal name starting with 'VCC' or 'VDD' is considered power.
* Any signal name starting with 'VSS' or 'GND' is considered ground.
* Added tests for various Inkscape path weirdness:
    * floating-point inaccuracy from relative moves
    * zero-length segments
    * shells that are clockwise (right-hand rule -> negative area)
    * self-crossing paths
    * multi-shell paths